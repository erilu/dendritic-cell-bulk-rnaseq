# Complete RNA seq pipeline: Transcriptome analysis of activated dendritic cells
Here, I will use **STAR** to align raw sequencing reads (.fastq files) to a reference genome, then use **Samtools** to count reads per gene. Last, use the **Bioconductor** package in **R** identify upregulated genes in activated dendritic cells.

---
## The Data

The raw data is publicly available for download at: https://www.ncbi.nlm.nih.gov/sra?term=SRP061381.
The GEO entry with details of the RNAseq experiment can be found at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1828772. There are four biological samples:
* two biological replicates of a sorted population of un-activated (control) dendritic cells
* two biological replicates of a sorted population of activated dendritic cells.

For each biological replicate, there will be several files. For this RNAseq experiment, paired-end sequencing was performed. The corresponding paired-end files are indicated by R1 and R2. If downloaded successfully, the raw data files should look like this when unzipped:

```
ls *.fastq.gz

# Het_CD4_1_AAGGGA_L001_R1_001.fastq.gz  IMM_CD4_1_TTCAGC_L001_R1_001.fastq.gz
# Het_CD4_1_AAGGGA_L001_R1_002.fastq.gz  IMM_CD4_1_TTCAGC_L001_R1_002.fastq.gz
# ...
# Het_CD4_1_AAGGGA_L001_R2_001.fastq.gz  IMM_CD4_1_TTCAGC_L001_R2_001.fastq.gz
# Het_CD4_1_AAGGGA_L001_R2_002.fastq.gz  IMM_CD4_1_TTCAGC_L001_R2_002.fastq.gz
# ...
# Het_CD4_2_AAGACG_L002_R1_001.fastq.gz  IMM_CD4_2_TTCGCT_L002_R1_001.fastq.gz
# ...
# Het_CD4_2_AAGACG_L002_R2_001.fastq.gz  IMM_CD4_2_TTCGCT_L002_R2_001.fastq.gz
# ...
```

Once the data is properly downloaded, we can begin the RNAseq analysis.
## The objective

When performing research in biology, one often needs to examine changes in gene expression in order to better understand why or how certain processes occur. RNAseq is very useful for determining the differences in global gene expression between two conditions, and is therefore particularly useful as a tool for discovery.

Our objective is to identify differentially expressed genes between un-activated and activated CD4+ dendritic cell populations. CD4+ dendritic cells are important for presenting antigen and activating CD4+ T cells. RNAseq can be used to identify markers and pathways that the dendritic cells might be using perform these functions.

---
# Aligning raw .fastq files to reference genome and counting reads

This data was previously processed using TopHat. An Excel file containing the differentially expressed genes is available as a supplementary on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71165).

We will process it using a different approach here, with the RNAseq alignment tool: [STAR](https://github.com/alexdobin/STAR). The code, along with an explanation of what each line does, is named [align_reads_STAR.sh](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/align_reads_STAR.sh).

---
# Performing differential gene analysis using DESeq2 (Bioconductor)

We can now read in the cleaned up data file and explore the differences in transcription between cell lines of interest. We will cluster the cell lines and pull out differentially expressed genes between subtypes of cells. To do this, we will use the [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/) package [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). This package has the capability to read in raw data files from each individual biological replicate, as well as the ability to analyze the read counts from an already compiled data file (which is what the Human Protein Atlas provides).

The file [cellatlas_analysis.R](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/cellatlas_analysis.R) will take you through how to perform the analysis (reading the data in, making a DESeq object, and annotating, exporting, and plotting the results). The exported list of differentially expressed genes between hematopoietic cells and non-hematopoietic cels is called [results_hemato_vs_non_DEGs.csv](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_hemato_vs_non_DEGs.csv). Below are some sample plots that can be used to visualize the data.

---
# Visualizations

### PCA clustering of cell lines
![PCA clustering of DC samples](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/results_PCA_plot.png)

Above is a principal components analysis (PCA) plot on a subset of the cell lines in the dataset. We can observe that the hematopoietic cell lines cluster away from the non-hematopoietic cell lines in the PCA plot. This suggests that they have different gene expression profiles. The next plot will show us some of the genes that contribute to the clustering we see here.

### Volcano plot to visualize differentially expressed genes with p-value cutoff

![Volcano plot cell line DEGs](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_volcano_plot_DEGs.png)

In accord with the clustering analysis, there are a lot of genes that are differentially expressed in hematopoietic vs non-hematopoietic cells.

### Heatmap to display top differentially expressed genes

![Heatmap of top DEGs](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis/blob/master/results_heatmap_top50_DEGs_ggplot2.png)

The differentially expressed genes can also be visualized using a heatmap.

---

If you are interested in learning how to perform a full RNA-seq pipeline analysis, you can look at my other [repo](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis) where I align raw .fastq sequencing files to a mouse reference genome, then use Bioconductor to find differentially expressed genes in activated vs. un-activated dendritic cells.
