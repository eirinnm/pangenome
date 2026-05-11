#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

counts_file <- args[1]
samplesheet <- read.csv(args[2])
out_file <- args[3]
ref_name <- args[4]

library(DESeq2)

counts <- read.delim(counts_file, comment.char="#", row.names=1)


# featureCounts adds annotation columns (Chr, Start, End, Strand, Length)
# Column 6 onwards are the BAM file paths
counts_data <- counts[, -(1:5), drop=FALSE]

# Get the current column names (these are file paths like 'work/xx/xxxx/A.bam')
current_bam_paths <- colnames(counts_data)

# Create a mapping by extracting the sample name from the file path
# This assumes your BAM files are named '${sample}.sorted.bam' as in your ALIGN process
extracted_samples <- gsub(paste0("_", ref_name, "\\.sorted\\.bam$"), "", basename(current_bam_paths))
# Reorder or verify the samplesheet matches the columns in the count matrix
# It's safest to subset the samplesheet to match the BAMs present in the matrix
samplesheet <- samplesheet[match(extracted_samples, samplesheet$sample), ]

# Now assign the clean sample names to the count matrix
colnames(counts_data) <- samplesheet$sample

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = samplesheet,
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), out_file)