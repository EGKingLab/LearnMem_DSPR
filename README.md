# Analyses from: Multiple genetic loci affect place learning and memory performance in Drosophila melanogaster

### Functions

The functions folder holds some general functions used in multiple scripts. 

1. mappingfunctions.R holds some basic functions used for QTL mapping
2. ggplot_theme holds the plotting theme for publication format, color scheme etc.
3. h2_functs.R holds functions for calculating heritabilities and genetic correlations

## Processing phenotypic data

Phenotype processing steps are in the LearnMemTolPheno folder. 

1. Process raw data for learning, memory, and thermal tolerance
  - Learning and Memory
      - Processed data from HeatCalc is parsed using a perl script (heatcalcRead_LearnMem.pl) and outputs: HeatProc_Learn.txt
      - Raw .asc files are processed using an R script (process_raw_learnmem.R) to identify any other   issues (see Step 2) and outputs Learn_raw.rda

2. Raw data and processed data are checked using check_raw_data.Rmd. Notes are within this script. Outputs: LearnMem_processed.txt

3. Processing_Pheno_data calculates the line means and outputs: L_MDATA.rda

4. Visualize_pheno_data makes the phenotype plots (Figure 2)

## QTL mapping and Quantitative Genetics 

QTL mapping analyses are in LearnMemTolQTL

1. mapping_perms.R performs the QTL mapping and permutations and calculation of FDR and FWER

2. identify_peaks.R finds each QTL peak and saves it in a list.

3. genes_under_peaks.R  looks at the DE genes from RNAseq (see below) to see which DE genes are under each peak

4. visualize_genome_scan.Rmd plots the genome scan (Figure 3)

5. QTL_effects_visualization plots the haplotype effects for each shared QTL peak (Figure 4)

6. IndividualPeaks.Rmd plots each peak interval separately 

7. Visualize_Peak_Intervals.Rmd plots each peak interval with the DE genes (Figure 6)

8. h2_rg.Rmd calculates heritabilities and genetic correlations


## RNA-Seq data analysis

### The alignment and assembly is in LearnMemTolRnaseq/RNAseq_Processing_LearnMemTol

1.Set_up_commands_masterpipeline.R is the master script that outputs 8 separate scripts for processing RNAseq data. It sets up scripts to run the RNA-Seq data on a cluster using SLURM. It will require adaptation for your own system and file paths.

Output sripts are named: 
S01_Trim_LearnMemRNA.txt - trims data
S02_QC_LearnMemRNA.txt - performs quality control
S03_Align_LearnMemRNA.txt - aligns data to the fly genome 
S04_SamtoBam_LearnMemRNA.txt - converts file from .sam to .bam
SO5_MergeBam_LearnMemRNA.txt - merges bam files
SO6_Assemble_LearnmemRNA.txt - assembles data

S08_Abundances_LearnMemRNA_eb.txt - finds abundances (Stringtie shorter version. This step only identifies genes that are already known. No novel genes are detected.)
--or--
S08_Abundances_LearnMemRNA.txt - finds abundances including novel transcripts


### Differential Expression Analysis is in LearnMemTolRnaseq/DEseq_LearnMemTol

1. prepDEseq.py takes the output form S08_Abundances_LearnMemRNA_pl.txt, stringtie (shorter version) and converts it into a usual version for the DEseq pipeline.
2. LearnMem_RNAproc_DESeq_eb.Rmd runs the DEseq analysis pipeline on the learning and memory RNAseq data.
3. Visualize_LMRNA_Data.Rmd plot graphs (Figure 5).

