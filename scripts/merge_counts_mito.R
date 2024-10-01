##
## Gazal Kalyan
## Jason Podrabsky Lab, Department of Biology
## Portland State Univerisity
## gazal@pdx.edu, podrabsj@pdx.edu
##

## To start with no saved variables and functions
## Uncomment next line
#rm(list = ls())

#####################################################
## Functions
#####################################################

# Change the following function  to get sample names from filenames
# These name will form column names  and
# Column names should match Row names from Sample Design table

samplename_from_filename <- function(filename){
    samplename <- unlist(strsplit(basename(filename), split = "\\."))[1]
    return(samplename)
}

get_input_filepath <- function(filepattern) {
    # Read input files
    input_files <- list.files(pattern = filepattern, 
                        recursive = TRUE)
    # Warnings
    if (length(input_files) == 0) {
        warning("No files fount for given filepattern")
    } else if (length(input_files) != NROW(samples)) {
        warning("Number of samples does not match the number of input files")
    }  else {
        warning("Successfully retrieved filepaths") }
    return(input_files)
}

make_filename_df <- function(filepaths, filetype) {
    for (i in 1:length(filepaths)) {
        samplename <- samplename_from_filename(filepaths[i])
        df_files_tmp <- data.frame(Sample = samplename)
        df_files_tmp[ ,filetype] = filepaths[i]
        if (i == 1) {
            df_files <- df_files_tmp
        } else {
            df_files <- rbind(df_files, df_files_tmp)
        }
        rm (df_files_tmp)
    }
    return(df_files)
}
#####################################################
## Set File Directory and File Paths
#####################################################
# The parent file directory that contains all result text files

parent_dir <- "/home/gazal/Documents/PSU_work/sRNA_gazal/OUTPUTS_FINAL/alim_mito_anno"
samplespath <- "samples.csv"

# Set pattern of filenames,
# these files are searched in the parent_dir recursively
# Make sure that the parent_dir contains these files and
# the pattern should be unique to list the files that we need


# 1. From **sports** annotation: *output.txt
# 2. From .bam files with **samtools**: *summary.txt
# 3. From .bam files with **samtools** and **awk**: *readcount.txt
# 4. From .bam files with **bedtools**: *bedtools.txt

txt_pattern_readcount <- "sorted.readcount.txt$"
txt_pattern_output <- "mapped_output.txt$"
txt_pattern_summary <- "txt.summary.txt$"
txt_pattern_bedtools <- "bedtools.txt$"

#####################################################
## Create Pheno table
#####################################################

library(dplyr)
library(data.table)

setwd(file.path(parent_dir))
samples <- read.csv(samplespath, sep = ",", header = TRUE)
rownames(samples) <- samples$Sample

print("Number of Samples found in the design file:")
print(NROW(samples))
print("Number of Observations found in the design file:")
print(NCOL(samples))
rm(samplespath)

#####################################################
## Create FilePath table
#####################################################

files_readcount <- get_input_filepath(txt_pattern_readcount)
files_output <- get_input_filepath(txt_pattern_output)
files_summary <- get_input_filepath(txt_pattern_summary)
files_bedtools <- get_input_filepath(txt_pattern_bedtools)

rm(txt_pattern_readcount,
    txt_pattern_output,
    txt_pattern_summary,
    txt_pattern_bedtools)

#Debug
############
print("Checking filenames")
sample_test <- samplename_from_filename(files_readcount[1])
print("The column names for each sample will be of the format:")
print(sample_test)
rm(sample_test)
#############

df_files_readcount <- make_filename_df(files_readcount, "readcountfile")
df_files_output <- make_filename_df(files_output, "outputfile")
df_files_summary <- make_filename_df(files_summary, "summaryfile")
df_files_bedtools <- make_filename_df(files_bedtools, "bedtoolsfile")

rm(files_readcount,
    files_output,
    files_summary,
    files_bedtools)

df_files <- merge(df_files_readcount, df_files_output, by = "Sample",
                all.x = TRUE, all.y = TRUE)
df_files <- merge(df_files, df_files_summary, by = "Sample",
                all.x = TRUE, all.y = TRUE)
df_files <- merge(df_files, df_files_bedtools, by = "Sample",
                all.x = TRUE, all.y = TRUE)

rm(df_files_readcount,
    df_files_output,
    df_files_summary,
    df_files_bedtools)

print(df_files)

#####################################################
## Create Readcount table from Bowtie alignments
#####################################################

for (i in 1:NROW(df_files)) {
    df_count_tmp <- read.csv(df_files$readcountfile[i],
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE,
        na.strings = ""
    )
    colnames(df_count_tmp) <- c("Counts", "Seq")
    df_count_tmp <- setorder(df_count_tmp, Seq)
    df_sum_tmp <- read.csv(df_files$summaryfile[i],
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE)
    colnames(df_sum_tmp) <- c("ReadID", "Seq", "Reference", "Start", "Alignments")
    df_sum_tmp <- df_sum_tmp[c("Seq", "Alignments")]
    df_sum_tmp$Alignments <- gsub("XM:i:", "", df_sum_tmp$Alignments)
    df_sum_tmp <- setorder(df_sum_tmp, Seq)
    df_sum_tmp <- distinct(df_sum_tmp, Seq, .keep_all = TRUE)

    df_sum_count_tmp <- merge(df_sum_tmp, df_count_tmp, by = "Seq", all.x = TRUE, all.y = TRUE)
    df_sum_count_tmp$AdjCount <- as.numeric(df_sum_count_tmp$Counts)/as.numeric(df_sum_count_tmp$Alignments)
    df_sum_count_tmp <- df_sum_count_tmp[df_sum_count_tmp$AdjCount > 1, ]
    df_sum_count_tmp <- df_sum_count_tmp[c("Seq", "AdjCount")]
    sample <- df_files$Sample[i]
    colnames(df_sum_count_tmp) <- c("Seq", sample)

    if (i == 1) {
        df_counts <- df_sum_count_tmp
    } else {
        df_counts <- merge(df_counts, df_sum_count_tmp, by = "Seq", all.x = TRUE, all.y = TRUE)
    }
    rm(df_sum_count_tmp, sample, df_count_tmp, df_sum_tmp)
    cat("readcount file:", i, " ")
}
rm(i)
df_counts <- distinct(df_counts, Seq, .keep_all = TRUE)
df_counts[is.na(df_counts)] <- 0
row.names(df_counts) <- df_counts$Seq
df_counts_new <- df_counts[, -1]
rm(df_counts)

# Remove low count sequences
df_counts_new <- df_counts_new[rowMeans(df_counts_new) > 2, ]
write.csv(df_counts_new, "out_tmp/output_readcount_table.csv")
df_counts_new <- read.csv("out_tmp/output_readcount_table.csv",
        sep = ",", row.names = 1, header = TRUE)
# TODO Divide by num of alignments

#####################################################
## Create Annotation table from Sports
#####################################################

for (i in 1:NROW(df_files)) {
    df_anno_tmp <- read.csv(df_files$outputfile[i],
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
    )
    colnames(df_anno_tmp)
    # Columns are ID, Sequence, Length, Reads, Match_Genome, Annotation
    # df_anno_tmp <- df_anno_tmp[df_anno_tmp$Reads > 0, ]
    df_anno_tmp <- df_anno_tmp[c("ID", "Sequence", "Length", "Annotation")]
    if (i == 1) {
        df_anno <- df_anno_tmp
    } else {
        df_anno <- rbind(df_anno, df_anno_tmp)
    }
    df_anno <- setorder(df_anno, Sequence)
    df_anno <- distinct(df_anno, Sequence, .keep_all = TRUE)
    rm(df_anno_tmp)
    cat("output file:", i, " ")
}
rm(i)
write.csv(df_anno, "out_tmp/output_annotation_table.csv", row.names = FALSE)
df_anno <- read.csv("out_tmp/output_annotation_table.csv",
        sep = ",", header = TRUE)

# TODO replace "Reads" column name to its sample name
# "Reads" column is redundant

#####################################################
## Create Summary table from Bowtie alignments
#####################################################
for (i in  1:NROW(df_files)) {
    bed <- read.csv(df_files$bedtoolsfile[i], 
        header = FALSE, 
        sep = "", 
        stringsAsFactors = FALSE)
    df_bed_tmp <- tidyr::separate(
        data = bed, col = V2, sep = ";",
        into = c("geneID", "Name", "gbkey", "gene", "gene_biotype", "gene_synonym")
    )
    rm(bed)
    df_bed_tmp <- df_bed_tmp[c("V1", "geneID", "gene", "gbkey", "gene_biotype")]
    df_bed_tmp$geneID <- gsub("ID=", "", df_bed_tmp$geneID)
    df_bed_tmp$gene <- gsub("gene=", "", df_bed_tmp$gene)
    df_bed_tmp$gbkey <- gsub("gbkey=", "", df_bed_tmp$gbkey)
    df_bed_tmp$gene_biotype <- gsub("gene_biotype=", "", df_bed_tmp$gene_biotype)
    colnames(df_bed_tmp) <- c("ReadID", "geneID", "gene", "gbkey", "gene_biotype")

    df_sum_tmp <- read.csv(df_files$summaryfile[i],
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE)
    colnames(df_sum_tmp) <- c("ReadID", "Seq", "Reference", "Start", "Alignments")
    df_sum_tmp <- df_sum_tmp[c("ReadID", "Seq", "Reference", "Start")]
    df_sum_tmp <- setorder(df_sum_tmp, ReadID)
    df_sum_tmp <- distinct(df_sum_tmp, ReadID, .keep_all = TRUE)

    df_bed_sum_tmp <- merge(df_sum_tmp, df_bed_tmp,
                    by = "ReadID", 
                    all.x=TRUE, 
                    all.y=FALSE)
    df_bed_sum_tmp <- setorder(df_bed_sum_tmp, Seq)
    df_bed_sum_tmp <- distinct(df_bed_sum_tmp, Seq, .keep_all = TRUE)
    rm(df_bed_tmp, df_sum_tmp)
    if (i == 1) {
        df_bed_sum <- df_bed_sum_tmp
    } else {
        df_bed_sum <- rbind(df_bed_sum, df_bed_sum_tmp)
    }
    rm(df_bed_sum_tmp)
    cat("summary file:", i, " ")
    df_bed_sum <- setorder(df_bed_sum, Seq)
    df_bed_sum <- distinct(df_bed_sum, Seq, .keep_all = TRUE)
}
rm(i)
write.csv(df_bed_sum, "out_tmp/output_bedtools__summary_table.csv", row.names = FALSE)
df_bed_sum <- read.csv("out_tmp/output_bedtools__summary_table.csv",
       sep = ",", header = TRUE)

#####################################################
## Create Location + Alignment table
#####################################################
# colnames are Seq
colnames(df_counts_new)
# Sequence column is unique
colnames(df_anno)
# ReadID is unique
colnames(df_bed_sum)

df_bed_sum_anno <- merge(df_anno, df_bed_sum,
    by.x = "Sequence", by.y = "Seq",
    all.x = TRUE, all.y = TRUE
)
df_bed_sum_anno <- setorder(df_bed_sum_anno, Sequence)
df_bed_sum_anno <- distinct(df_bed_sum_anno, Sequence, .keep_all = TRUE)
row.names(df_counts) <- df_counts$Sequence
write.csv(df_bed_sum_anno, "out_tmp/output_anno_bed_sum_table.csv")

df_bed_sum_anno_count <- merge(df_counts_new, df_bed_sum_anno,
    by.x = "row.names", by.y = "Sequence",
    all.x = TRUE, all.y = FALSE
)
write.csv(df_bed_sum_anno_count, "out_tmp/output_anno_bed_sum_counts_table.csv")
