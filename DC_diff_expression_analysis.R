
# DC_diff_expression_analysis.R
# After running STAR pipeline on raw fastq files, the output was an Htseq file per biological sample
# This code will analyze differential gene expression between unimmunized and immunized dendritic cells (DCs)
# using the bioconductor package DESeq2, and visualize the data using ggplot2

# Initialize the packages required to analyze RNAseq data
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2") 
library("ggrepel") 
library("gplots")

# set the working directory into the directory that contains your htseq files.
htseq_dir = "~/Bioinformatics/DCrnaseq"
setwd(htseq_dir)

# find all the htseq files in the directory and list them, naming them by what they represent (WT or IMM)
htseq_files = list.files(pattern = "*_htseq.out")
names(htseq_files) = c ("Control_1", "Control_2", "Activated_1", "Activated_2")

# Table should look like this:
# Control_1                    Control_2                    Activated_1                  Activated_2 
# "Het_CD4_1_byName_htseq.out" "Het_CD4_2_byName_htseq.out" "IMM_CD4_1_byName_htseq.out" "IMM_CD4_2_byName_htseq.out"

# Read in the data from the the htseq files into a table format. This is the format required for the DESeq function, countDataSet
# If you have many different samples, gsub() can pull out hte prefix
condition = factor(gsub("_\\d+$", "", names(htseq_files)) )
sampleTable = data.frame(sample_name = names(htseq_files), file_name = htseq_files, condition = relevel(condition, ref = "Control"))

# Sample table should look like this:
#             sample_name                  file_name condition
# Control_1     Control_1 Het_CD4_1_byName_htseq.out   Control
# Control_2     Control_2 Het_CD4_2_byName_htseq.out   Control
# Activated_1 Activated_1 IMM_CD4_1_byName_htseq.out Activated
# Activated_2 Activated_2 IMM_CD4_2_byName_htseq.out Activated

# Run the DESeq package's function, newCountDataSetFromHTSeqCount
# 
dds <- DESeqDataSetFromHTSeqCount(sampleTable, directory = ".", design = ~ condition)

# can perform regularized-log transformation before plotting data on PCA plot:
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)
save(rld, file = "rld_dds.Robj")

# Use the rlog transformed data to cluster the cell lines based on similarity
plotDists = function (rld.obj) {
  sampleDists <- dist(t(assay(rld.obj)))
  
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste( rld.obj$celltype )
  
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
}

plotDists(rld)

# https://support.bioconductor.org/p/90791/
# plotPCA() shows the PCA plot of the log transformed data
# This code adds the names of the samples ontop of the PCA plot dots
name.plotPCA = function (rld.obj) {
  p <- plotPCA(rld.obj,  intgroup = c("condition"))
  p <- p + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
  p = p+ ggtitle("PCA Plot - Unactivated vs 1 hr SRBC-activated dendritic cells")
  print(p)
}

name.plotPCA(rld)

# Run DESeq() on the dataset to test for differentially expressed genes between the two conditions

dds_deseq = DESeq(dds)
res = results(dds_deseq)

# This is what the res table should look like at this step:

# log2 fold change (MAP): condition Activated vs Control 
# Wald test p-value: condition Activated vs Control 
# DataFrame with 6 rows and 6 columns
#                 baseMean  log2FoldChange     lfcSE        stat       pvalue         padj
#                 <numeric>      <numeric> <numeric>   <numeric>    <numeric>    <numeric>
# 0610005C13Rik   3.5593163    -0.59338435 0.5564591 -1.06635755 2.862620e-01           NA
# 0610007C21Rik 151.3048006    -0.76444801 0.2804598 -2.72569571 6.416612e-03 2.885478e-02
# 0610007L01Rik 542.9491558     0.83896607 0.1796740  4.66938003 3.021101e-06 3.349804e-05
# 0610007N19Rik   0.8000586     0.07346235 0.3350871  0.21923358 8.264681e-01           NA
# 0610007P08Rik 208.1359047    -0.50145544 0.2479413 -2.02247632 4.312717e-02 1.290959e-01
# 0610007P14Rik 220.4997078    -0.01215840 0.2394232 -0.05078207 9.594992e-01 9.780066e-01

# append the raw counts to the results table so you can scrutinize the data per biological replicate
raw_counts = assay(dds)
resTable = cbind(res, raw_counts)

# we want to annotate the files with full gene games and Gene Ontology terms, using the mouse genes annotation package
library("AnnotationDbi")
library("org.Mm.eg.db")

# This function adds columns representing gene ontology term and the expanded gene name
my.mapids = function (res) {
  res$go <- mapIds(org.Mm.eg.db,
                   keys=row.names(res),
                   column="GO",
                   keytype="SYMBOL",
                   multiVals="first")
  res$genename <- mapIds(org.Mm.eg.db,
                         keys=row.names(res),
                         column="GENENAME",
                         keytype="SYMBOL",
                         multiVals="first")
  resOrdered <- res[order(res$padj),]
  return(resOrdered)
}

res_annotated = my.mapids(resTable)
# save the RObjects to files incase we want to load them for analysis later on.
save(res_annotated, file = "res_annotated.Robj")

# export gene lists, top 2000 differentially expressed genes, sorted by padj
resOrderedDF <- as.data.frame(res_annotated)[1:2000, ]

write.csv(resOrderedDF, file = "top2000_diff_genes_padj.csv")

# alternatively we could cut off all genes below a certain p value and return an ordered list of the biggest differences

sort_diff_genes = function(df, numberofgenes=200) {
  df = df[which(df$padj < 0.0001),]
  df = df[order(abs(df$log2FoldChange), decreasing=TRUE)[1:numberofgenes],]
  df = df[order(df$log2FoldChange, decreasing=TRUE ), ]
  return(df)
}

sorted_res = sort_diff_genes(res_annotated,200)
write.csv(sorted_res, file = "top200_diff_genes_log2fchange.csv")


# Volcano plot to visualize top 20 differentially expressed genes in the cleaned up dataset

plot.volcano = function (res) {
  input <- mutate(data.frame(res), significance=ifelse(data.frame(res)$padj<0.0001, "padj < 0.0001", "Not Significant"))
  
  input = input[!is.na(input$significance),]
  volc = ggplot(input, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(col=significance)) + 
    scale_color_manual(values=c("black", "red")) + 
    ggtitle("Volcano Plot - Differentially Expressed Genes 
            \nEnriched in activated (right) vs unactivated (left) CD4 dendritic cells") +
    geom_text_repel(data=head(input, 20), aes(label=rownames(head(res,20))))
  print(volc)
}

plot.volcano(res_annotated)

# plot a heatmap using ggplot2

# ggplot2 heatmap requires long format data - we currently have data in the "wide" format -
# to reorganize data, use melt - https://stackoverflow.com/questions/30040420/heat-map-per-column-with-ggplot2
# long vs wide data format - https://www.theanalysisfactor.com/wide-and-long-data/ 

rawmatrix = data.matrix(scale(sorted_res[c(1:20),-c(1:6,11,12)]))
melted = melt(rawmatrix)
colnames(melted) = c("Gene","Sample","Expression")

# scale each row individually (can also use dplyr, or apply() on the rawmatrix:
# https://stackoverflow.com/questions/33642723/scale-rows-of-data
# scale1 = apply(rawmatrix, 1, scale)
# melt_ver1 = melt(scale1)
# colnames(melt_ver1) = c("Gene","Sample","Expression")


scale_melted <- ddply(melted, .(Gene), transform, rescale = scale(Expression))

# The order of the differentially expressed genes matters in this case (highest log2foldchange)
# We want to preserve the order for the graph, so we will define a fixed order to them.
# https://rstudio-pubs-static.s3.amazonaws.com/7433_4537ea5073dc4162950abb715f513469.html
# https://stackoverflow.com/questions/8713462/ggplot2-change-order-of-display-of-a-factor-variable-on-an-axis
# https://stackoverflow.com/questions/36651511/order-data-for-ggplot2-geom-tile


scale_melted$Gene <- factor(scale_melted$Gene, levels = unique(rev(as.character(scale_melted$Gene))))

heatmap = ggplot(scale_melted, aes (Sample, Gene) ) +
  geom_tile(aes(fill = rescale), color = "white") +
  scale_fill_gradient(low = "white", high = "firebrick4") +
  ylab("List of Genes") +
  xlab("Samples") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  labs(fill = "Expression level") +
  ggtitle("Heatmap- Upregulated genes in \nunactivated vs 1 hr SRBC-activated dendritic cells")

print(heatmap)



# alternatively, can plot a heatmap using heatmap.2 from gplots package. formatting is difficult so
# you could export the plot as a PDF and then format with adobe illustrator. ggplot2 is more customizable

# this heatmap.2 requires data.matrix object, initialize one that has hematopoietic cell lines grouped together
raw = data.matrix(sorted_res[c(1:20),-c(1:6,11,12)])

# initialize variable for RColorBrewer (value can be either: BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral)
brewer_palette <- "RdBu"

# Ramp the color in order to get the scale.
ramp <- colorRampPalette(brewer.pal(11, brewer_palette))
mr <- ramp(256)[256:1]

# Run the heatmap.2 function in order to generate the heatmap.
# scale = "row" means each row will be standardized, "none" means raw values used
# col = mr uses the color scheme that we ramped from colorbrewer
# rowv / colv = FALSE means that the original order of the rows / cols will be retained.
# lwid and lhei change the relative size of the heatmap components

heatmap.2(raw,
          density.info="none",
          trace="none",
          margins =c(5,5),
          col=mr,
          scale = "row",
          dendrogram = 'none',
          srtCol = 45,
          cexCol = 1,
          lwid = c(0.05,4),
          lhei = c(1,5),
          main = 'Heatmap of top 20 upregulated genes between\nactivated and unactivated CD4 dendritic cells'
) 



