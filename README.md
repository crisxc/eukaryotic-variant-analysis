# Eukaryotic Variant Analysis

Scripts for whole-genome sequencing variant calling, filtering, and mutation analysis.
This has been tested on fungal Candida species, but can be adapted for any other organism with a reference genome.

## Workflow

1. Call variants
    
    Generates colony-level variant calls from paired FASTQs using GATK joint genotyping.
    Shoutout to Elizabeth Mei for integrating GATK steps!   
2. Analyze variants
    
    Processes called variants for downstream analysis, including background strain filtering, colony-unique variant identification, & classification by mutation type.
    
3. Create figures
    
    Visualize mutation burden, distribution, and chromosome-level patterns. Likely to be updated with more figure types in the future.
    

## Folders

### 01_call_variants/

Variant calling workflow starting from FASTQ files.

Includes:

- read alignment (FASTQ input)
- GATK Haplotypecaller
- joint genotyping across colony samples
- variant filtering
- TSV export for step 2

### 02_analyze_variants/

Processes joint-genotyped variant tables after calling.

Includes:

- background strain variant filtering
- colony-unique mutation identification
- mutation type classification
- strain-, colony-, and chromosme-level mutation count tables

### 03_create_figures/

Scripts for visualizing mutation patterns.

Includes:

- chromosome mutation maps
- mutation burden plots
- genomic distribution plots
