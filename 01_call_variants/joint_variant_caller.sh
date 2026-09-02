#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

usage() {
    echo "Usage: $0 -f FASTQ_DIR -r reference.fasta -o output_dir [-t threads] [-n run_name]"
    exit 1
}

fastq=""
ref=""
out="joint_variant_results"
threads=8
run_name=""

while getopts "f:r:o:t:n:" opt; do
    case "$opt" in
        f) fastq="$OPTARG" ;;
        r) ref="$OPTARG" ;;
        o) out="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        n) run_name="$OPTARG" ;;
        *) usage ;;
    esac
done

[[ -n "$fastq" ]] || usage
[[ -n "$ref" ]] || usage

fastq=$(realpath "$fastq")
ref=$(realpath "$ref")
out_parent=$(realpath -m "$out")

first_r1=$(find "$fastq" -maxdepth 1 -type f -name "*R1*fastq*" | sort | head -n1)
[[ -n "$first_r1" ]] || { echo "Error!! no R1 fastq found in $fastq"; exit 1; }

if [[ -z "$run_name" ]]; then
    run_name=$(basename "$first_r1" | sed -E 's/(.*)[-_]R1[-_]?.*/\1/')
fi

mkdir -p "$out_parent/${run_name}_variant_calls"
out=$(realpath "$out_parent/${run_name}_variant_calls")

prefix="$out/$run_name"
log="${prefix}_pipeline.log"

{
    echo "Started: $(date)"
    echo "FASTQ: $fastq"
    echo "Reference: $ref"
    echo "Output: $out"
    echo "Run name: $run_name"
} > "$log"

ref_base=$(basename "$ref")
ref_base="${ref_base%.gz}"
ref_base="${ref_base%.fasta}"
ref_base="${ref_base%.fa}"

ref_dir=$(dirname "$ref")
ref_dict="$ref_dir/${ref_base}.dict"
idx="$ref_dir/$ref_base"

# reference files
[[ -f "$ref.fai" ]] || samtools faidx "$ref" >> "$log" 2>&1
[[ -f "$ref_dict" ]] || gatk CreateSequenceDictionary -R "$ref" -O "$ref_dict" >> "$log" 2>&1

# bowtie2 index
if [[ ! -f "$idx.1.bt2" && ! -f "$idx.1.bt2l" ]]; then
    bowtie2-build "$ref" "$idx" >> "$log" 2>&1
fi

cd "$out"

# FASTQ to BAM
for r1 in "$fastq"/*R1*fastq*; do
    base=$(basename "$r1")
    sample=$(echo "$base" | sed -E 's/(.*)[-_]R1[-_]?.*/\1/')
    r2=$(ls "$fastq/${sample}"*R2*fastq* 2>/dev/null | head -n1 || true)
    sample_log="$out/${sample}_sample.log"

    [[ -z "$r2" ]] && { echo "missing R2 for $sample" >> "$log"; continue; }
    [[ -f "${sample}_RG_MD.bam" ]] && { echo "$sample BAM exists; so i'm skipping it" >> "$log"; continue; }

    echo "aligning $sample" >> "$log"
    echo "Started: $(date)" > "$sample_log"

    bowtie2 -x "$idx" -1 "$r1" -2 "$r2" --threads "$threads" \
        -S "${sample}.sam" >> "$sample_log" 2>&1

    samtools view -bS "${sample}.sam" | \
        samtools sort -@ "$threads" -o "${sample}_sorted.bam" >> "$sample_log" 2>&1

    gatk AddOrReplaceReadGroups -I "${sample}_sorted.bam" -O "${sample}_RG.bam" \
        --RGID "$sample" --RGLB "$sample" --RGPL ILLUMINA --RGPU unit1 --RGSM "$sample" >> "$sample_log" 2>&1

    gatk MarkDuplicates -I "${sample}_RG.bam" -O "${sample}_RG_MD.bam" \
        --METRICS_FILE "${sample}_metrics.txt" >> "$sample_log" 2>&1

    samtools index "${sample}_RG_MD.bam" >> "$sample_log" 2>&1
    rm -f "${sample}.sam" "${sample}_sorted.bam" "${sample}_RG.bam"

    echo "Done: $(date)" >> "$sample_log"
done

# BAM to gVCF
for b in *_RG_MD.bam; do
    sample=$(basename "$b" _RG_MD.bam)
    sample_log="$out/${sample}_sample.log"

    [[ -f "${sample}.g.vcf.gz" ]] && { echo "$sample gVCF exists; skipping!" >> "$log"; continue; }

    echo "calling $sample" >> "$log"

    gatk HaplotypeCaller -R "$ref" -I "$b" \
        -O "${sample}.g.vcf.gz" -ERC GVCF \
        -G StandardAnnotation >> "$sample_log" 2>&1
done

# combine samples
gvcfs=""
for g in *.g.vcf.gz; do
    gvcfs="$gvcfs -V $g"
done

[[ -n "$gvcfs" ]] || { echo "Error!! no gVCFs were made. Check log file."; exit 1; }

gatk CombineGVCFs -R "$ref" $gvcfs \
    -O "${prefix}.g.vcf.gz" >> "$log" 2>&1

gatk GenotypeGVCFs -R "$ref" -V "${prefix}.g.vcf.gz" \
    -O "${prefix}_genotyped.vcf.gz" >> "$log" 2>&1

# filter variants. Change if necessary.
gatk VariantFiltration -R "$ref" -V "${prefix}_genotyped.vcf.gz" \
    -O "${prefix}_filtered.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name QD2 \
    --filter-expression "QUAL < 50.0" --filter-name QUAL50 \
    --filter-expression "FS > 60.0" --filter-name FS60 \
    --filter-expression "MQ < 40.0" --filter-name MQ40 \
    --filter-expression "MQRankSum < -12.5" --filter-name MQRankSum-12.5 \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRankSum-8 \
    --filter-expression "DP < 20.0" --filter-name DP20 >> "$log" 2>&1

gatk SelectVariants -R "$ref" -V "${prefix}_filtered.vcf.gz" \
    --exclude-filtered -O "${prefix}_filtered_PASS.vcf.gz" >> "$log" 2>&1

gatk SelectVariants -V "${prefix}_filtered_PASS.vcf.gz" -select-type SNP \
    -O "${prefix}_filtered_SNPs.vcf.gz" >> "$log" 2>&1

gatk SelectVariants -V "${prefix}_filtered_PASS.vcf.gz" -select-type INDEL \
    -O "${prefix}_filtered_INDELs.vcf.gz" >> "$log" 2>&1

# export table
gatk VariantsToTable -V "${prefix}_filtered_PASS.vcf.gz" \
    -F CHROM -F POS -F TYPE -F REF -F ALT -F AC -F AF \
    -F EVENTLENGTH -F TRANSITION -F HOM-REF -F HOM-VAR \
    -F HET -F VAR -F NO-CALL -F NCALLED -GF GT \
    -O "${prefix}_filtered_variants.tsv" >> "$log" 2>&1

echo "Done: $(date)" >> "$log"
