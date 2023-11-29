# Install and load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)

# Load count matrices
gene_counts <- read.csv("gene_count_matrix.csv", row.names = 1)
transcriptome_counts <- read.csv("transcript_count_matrix.csv", row.names = 1)

# Adding small count to address zero gene-counts issue
gene_counts <- gene_counts + 1
transcriptome_counts <- transcriptome_counts + 1

# Filter out low-count genes
gene_counts_filtered <- gene_counts[rowSums(gene_counts) > 1, ]
transcriptome_counts_filtered <- transcriptome_counts[rowSums(transcriptome_counts) > 1, ]

# Sample information (replace with your actual sample information)
col_data <- data.frame(
  cell_line = rep(c("U251", "U343"), each = 8),  # Assuming each cell line has 9 samples
  time = rep(c("T1", "T2"), times = 8)  # Assuming each time point has 9 samples
)

# print the deimensions of matrix
print(dim(gene_counts))
print(dim(transcriptome_counts))
print(dim(col_data))

# Create DESeqDataSet objects
dds_gene <- DESeqDataSetFromMatrix(countData = gene_counts_filtered, colData = col_data, design = ~ cell_line + time + cell_line:time)
dds_transcriptome <- DESeqDataSetFromMatrix(countData = transcriptome_counts_filtered, colData = col_data, design = ~ cell_line + time + cell_line:time)

# Run DESeq analysis
dds_gene <- DESeq(dds_gene)
dds_transcriptome <- DESeq(dds_transcriptome)

# Perform time-specific contrasts
contrast_gene <- c("time", "T2", "T1")
contrast_transcriptome <- c("time", "T2", "T1")

dds_gene_results <- results(dds_gene, contrast = contrast_gene)
dds_transcriptome_results <- results(dds_transcriptome, contrast = contrast_transcriptome)

# Adjust p-values
dds_gene_results <- results(dds_gene, contrast = contrast_gene, alpha = 0.05)
dds_transcriptome_results <- results(dds_transcriptome, contrast = contrast_transcriptome, alpha = 0.05)

# Filter differentially expressed genes
DE_genes <- subset(dds_gene_results, padj < 0.05)
DE_transcripts <- subset(dds_transcriptome_results, padj < 0.05)

# Explore and visualize results (example: MA plot)
plotMA(dds_gene_results, main="DESeq2", ylim=c(-2,2))

# Save results
write.csv(DE_genes, "differential_genes.csv")
write.csv(DE_transcripts, "differential_transcripts.csv")