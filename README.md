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

Our objective is to identify differentially expressed genes between un-activated and activated CD4+ dendritic cell populations. CD4+ dendritic cells are important for presenting antigen and activating CD4+ T cells. RNAseq can be used to identify markers and pathways that the dendritic cells might be using perform these functions. The results of this analysis containing a list of the top differentially expressed genes is here: [results_top2000_diff_genes_padj.csv]().

These **.fastq** files were previously processed using an alternate RNAseq analysis tool called TopHat. An Excel file containing the differentially expressed genes is available as a supplementary on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71165). The previously processed data was used in a [research study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4883664/) that I contributed to. The analysis that I walk through here is an updated version of the analysis used in the research paper using newer tools. Comparing the output from the two different methods, you can observe that both methods are identifying similar differentially expressed genes.

---
# Aligning raw .fastq files to reference genome and counting reads

We will process the data using the RNAseq alignment tool: [STAR](https://github.com/alexdobin/STAR). The code, along with an explanation of what each line does, is named [align_reads_STAR.sh](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/align_reads_STAR.sh).

The code first uses regexes to identify all the files that correspond to a particular biological replicate, then passes the files to the STAR function:

```
STAR --runThreadN 8 --genomeDir $genomedir/
 --genomeLoad LoadAndKeep --readFilesIn $end1 $end2
  --readFilesCommand zcat --outFileNamePrefix $prefix
   --outStd SAM --outFilterMultimapNmax 1 > ${prefix}.sam
```

This spits out a **.sam** file, which is sorted and then counted using ```htseq-count```. I've included the htseq-count output files in a [folder](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/tree/master/htseq_count_files) in this repo, in case you want to skip the alignment and analyze them in R yourself.


---
# Performing differential gene analysis using DESeq2 (Bioconductor)

Once we have the count files (\_htseq.out), we can read them in using R and perform differential gene analysis. [Bioconductor](https://www.bioconductor.org/help/workflows/rnaseqGene/) has a package called [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) that is commonly used for this purpose. This package has the capability to read in raw data files from each individual biological replicate (the \_htseq.out files), as well as analyze the read counts from an already compiled data file (if you get data from outside labs, it might be provided in this format). The file [DC_diff_expression_analysis.R](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/DC_diff_expression_analysis.R) contains the code used to perform the analysis.

The code will take you through the following analysis steps:
* reading the data in
* making a DESeq object
* creating a results object comparing two groups
* annotating, exporting, and plotting the results

The exported list of differentially expressed genes between unactivated and activated dendritic cells is called [results_top2000_diff_genes_padj.csv](). Below are some sample plots that can be used to visualize the data.

---
# Visualizations

### PCA clustering of cell lines
![PCA clustering of DC samples](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/results_PCA_plot.png)

Above is a principal components analysis (PCA) plot of the 4 biological samples. As expected, we can observe that the activated dendritic cell samples cluster away from the unactivated cell lines in the PCA plot. This suggests that the groups have different gene expression profiles. The next plot will show us some of the genes that contribute to the clustering we see here.

### Volcano plot to visualize differentially expressed genes with p-value cutoff

![Volcano plot of upregulated DEGs](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/results_volcano_plot_DEGs.png)

We can see that the majority of differentially expressed genes are upregulated in activated dendritic cells. This makes sense biologically, since a dendritic cell would want to rapidly begin expressing markers that help activate T cells (co-stimulatory molecules), such as CD86. CD86 is a well described marker that is upregulated by dendritic cells after activation. We can see that CD86 shows up in our analysis, which gives us more confidence that the experiment and data analysis was performed correctly.

### Heatmap to display top differentially expressed genes

![Heatmap of top DEGs](https://github.com/erilu/Complete-RNA-seq-Pipeline-Transcriptome-Analysis/blob/master/results_heatmap_DEGs.png)

The differentially expressed genes can also be visualized using a heatmap. I made an different version of this heatmap using the gene list from the previous analysis in [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71165), which was published in a [scientific research paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4883664/) in the journal _Nature_ (Extended Figure 5, panel A).

---

Thanks for reading! If you are interested in looking at another RNAseq analysis (using cell line RNAseq data obtained from the Human Protein Atlas), you can check out my other  [repo](https://github.com/erilu/R-Cell-Line-Transcriptome-Analysis). I clean and organize read count files, then use Bioconductor to find differentially expressed genes in subsets of cell lines. I also pull out all annotated enzymes and repeat the analysis to find differentially expressed enzymes.
