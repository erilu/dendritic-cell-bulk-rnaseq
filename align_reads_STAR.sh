# align_reads_STAR.sh
# This code reads in raw RNA-seq .fastq files, aligns them to a reference genome, and counts reads per gene using HTSeq.
# Raw data is publicly available for download at: https://www.ncbi.nlm.nih.gov/sra?term=SRP061381
# The GEO entry can be found at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1828772
# This data was previously analyzed using TopHat. We will re-analyze it here using STAR

###########################################################################################
# Setting up directories
###########################################################################################

# Download the raw fastq files from NCBI and store it in a directory of your choice.

# The raw data files should look like this (several files for each biological replicate, paired ends indicated by R1 and R2):
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

# Store the directory that contains the raw fastq files into a variable:
seqdir=/bigdata/seq/RnaSeq/

# Second, define the directory that will contain the analysis. If this doesn't exist, make it using "mkdir":
analysisdir=/bigdata/analysis/

# Third, define the directory that contains our reference genome (mm10).
genomedir=/bigdata/genomes/mm10/STAR

# Next, create soft-links to the raw data so that we can call the files in the analysis folder without typing out the whole directory each time.
# I like to separate the raw data from the analysis files so that the raw data is easier to locate in the future.
cd $analysisdir
find $seqdir -name '*.fastq.gz' | xargs -I {} sudo ln -s {} .

###########################################################################################
# Align your raw data to the reference genome using STAR
# Install STAR at: https://github.com/alexdobin/STAR
###########################################################################################

# The files have names like this: Het_CD4_1_AAGGGA_L001_R1_001.fastq.gz
# All the fastq files can be returned using this command: *.fastq.gz
# Here, "Het_CD4_1" is the name of one biological replicate. Each biological replicate has many files associated with it.
# Given these file names, we can use regexes to identify files that belong to a specific biological replicate.

mylist="$(ls *.fastq.gz | perl -pe 's[^([^_]+_[^_]+\d_\d)_[ACGT]{6}_L00\d_R\d_\d\d\d.fastq.gz][$1]' | uniq )"
echo $mylist

# This regex should return:
# Het_CD4_1 Het_CD4_2 IMM_CD4_1 IMM_CD4_2

# We can now identify all files that correspond to a biological replicate.
# For each replicate, we will merge/align all the corresponding files to a reference genome using the STAR command.
# The data I am analyzing is paired-end read data, so I must feed both paired ends into the STAR function at the same time.
# The paired ends are distinguished by the "R1" and "R2" in the file names:
# For example: Het_CD4_1_AAGGGA_L001_R1_001.fastq.gz and Het_CD4_1_AAGGGA_L001_R2_001.fastq.gz are paired end files.

# I like to print out what step the loop is currently on, to have an idea of when it will finish.
# I would also suggest testing the loop without running the actual STAR command, to ensure that you don't have a bug in your code.

for prefix in $mylist; do
end1=`ls ${prefix}_*_R1_*.fastq.gz | perl -e '@a=<>;chomp @a; print join q[,],@a'`
end2=`ls ${prefix}_*_R2_*.fastq.gz | perl -e '@a=<>;chomp @a; print join q[,],@a'`
echo ${prefix}
echo $end1
STAR --runThreadN 8 --genomeDir $genomedir/ --genomeLoad LoadAndKeep --readFilesIn $end1 $end2 --readFilesCommand zcat --outFileNamePrefix $prefix --outStd SAM --outFilterMultimapNmax 1 > ${prefix}.sam
echo 'aligned'
done

# What you should see:
# Het_CD4_1
# Het_CD4_1_AAGGGA_L001_R1_001.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_002.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_003.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_004.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_005.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_006.fastq.gz,Het_CD4_1_AAGGGA_L001_R1_007.fastq.gz
# aligned
# ...

# Running the STAR command generates the following: ${Prefix}Log.final.out, Log.out, Log.progress.out, Log.std.out, SJ.out.tab, and .sam.
# The file that contains the mapped reads is the .sam file.

# Note on single-end read data:
# Single-end read data is easier to work with, because each biological replicate is usually contained within its own fastq.gz file
# Since those files are not split up, you can just loop through each fastq file and directly align it to the genome using STAR.
# If the data presented here was single-end, below would be an example loop to align it. 
# You can remove the pipe "| perl -e '@a=<>;chomp @a; print join q[,],@a'`" if there is only one file per biological replicate.

# for prefix in $mylist; do
# end1=`ls ${prefix}_*_R1_*.fastq.gz | perl -e '@a=<>;chomp @a; print join q[,],@a'`
# echo $end1
# STAR --runThreadN 8 --genomeDir $genomedir/ --genomeLoad LoadAndKeep --readFilesIn $end1 --readFilesCommand zcat --outFileNamePrefix $prefix --outStd SAM --outFilterMultimapNmax 1 > ${prefix}.sam
# done

###########################################################################################
# Sorting the aligned output and processing raw read count per gene
###########################################################################################

# Now that you have your .sam files, you should sort reads by name for easier data analysis
# This function pipes the list of sam files ("ls *.sam") into samtools' "view" function, then pipes this samtools list into the "sort"
# function, which then sorts by readnames (-n), number of bytes of maximum memory (-m2GB), and specifies the output ext. ._byName
ls *.sam | parallel 'samtools view -bS {} | samtools sort -n -m2000000000 - {.}_byName'

# count mapped reads in each sam file and store this data in a text file
ls *.sam | parallel -k 'echo {}; grep -v ^@ {} | wc -l' > read_counts.txt

# make a soft link to the reference genome annotation file. This will be used by the next function.
ln -s /bigdata/genomes/mm10/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2012-05-23-16-47-35/Genes/genes.gtf .

# count the reads mapped to each gene (produces '_htseq.out' files)
# note: The default for strandedness is yes. If your RNA-Seq data has not been made with a strand-specific protocol, 
# this causes half of the reads to be lost. Hence, make sure to set the option --stranded=no unless you have strand-specific data!
# more info at: http://www-huber.embl.de/users/anders/HTSeq/doc/count.html
# What this function does is find all the _byName files (that you generated above), feed them into the samtools,
# which uses the htseq tool to count individual transcript reads for each of the genes in the file "genes.gtf"
find . -name '*_byName.bam' | parallel 'samtools view -h {} | htseq-count --stranded=no - genes.gtf > {.}_htseq.out'

# the resulting output files (*_htseq.out) are your raw sequence counts! You can now feed this file into R to
# perform differential gene expression analysis.
# for the code used to analyze these count files, refer to "DC_diff_expression_analysis.R"



