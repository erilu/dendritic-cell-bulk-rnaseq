
# DC_diff_expression_analysis.R
# After running STAR pipeline on raw fastq files, the output was an Htseq file per biological sample
# This code will analyze differential gene expression between unimmunized and immunized dendritic cells (DCs)
# using the bioconductor package DESeq2, as well as try computing the values manually using a negative bionmial model

# Initialize the packages required to analyze RNAseq data
library("DESeq2")
library("xtable")
library("parallel")
library("sqldf")

# set the working directory into the directory that contains your htseq files.
htseq_dir = "~/Desktop/Bioinformatics/DCrnaseq/htseq_count_files"
setwd(htseq_dir)

# find all the htseq files in the directory and list them, naming them by what they represent (WT or IMM)
htseq_files = list.files(pattern = "*_htseq.out")
names(htseq_files) = c ("Control_1", "Control_2", "Activated_1", "Activated_2")

# Table should look like this:
# Control_1                    Control_2                    Activated_1                  Activated_2 
# "Het_CD4_1_byName_htseq.out" "Het_CD4_2_byName_htseq.out" "IMM_CD4_1_byName_htseq.out" "IMM_CD4_2_byName_htseq.out"

# Read in the data from the the htseq files into a table format. This is the format required for the DESeq function, countDataSet
# If you have many different samples, gsub() can pull out hte prefix
condition = gsub("_\\d+$", "", names(htseq_files)) 
relevel (condition, ref = "Control")
sampleTable = data.frame(sample_name = names(htseq_files), file_name = htseq_files, condition = relevel (condition, ref = "Control"))

# Sample table should look like this:
#             sample_name                  file_name condition
# Control_1     Control_1 Het_CD4_1_byName_htseq.out   Control
# Control_2     Control_2 Het_CD4_2_byName_htseq.out   Control
# Activated_1 Activated_1 IMM_CD4_1_byName_htseq.out Activated
# Activated_2 Activated_2 IMM_CD4_2_byName_htseq.out Activated

# Run the DESeq package's function, newCountDataSetFromHTSeqCount
# 
dds <- DESeqDataSetFromHTSeqCount(sampleTable, directory = ".", design = ~ condition)

# Run DESeq() on the dataset to test for differentially expressed genes between the two conditions

dds_deseq = DESeq(dds)
res = results(dds_deseq)
resTable = cbind(res, assay(dds))


