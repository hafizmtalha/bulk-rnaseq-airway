# Import Salmon results into R and test for differential expression with DESeq2.
# The design accounts for the cell line each sample came from, so the model
# compares treated and control cells within the same line.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
})

samples <- read_tsv(snakemake@input[["samples"]], show_col_types = FALSE)
samples$condition <- factor(samples$condition, levels = c("control", "dex"))
samples$cell_line <- factor(samples$cell_line)

files <- file.path("results/salmon", samples$sample, "quant.sf")
names(files) <- samples$sample

tx2gene <- read_tsv(snakemake@input[["tx2gene"]], show_col_types = FALSE)
gene_info <- tx2gene %>% distinct(gene_id, gene_name)

# Sum transcript level estimates to genes
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("transcript_id", "gene_id")],
                dropInfReps = TRUE)

coldata <- as.data.frame(samples)
rownames(coldata) <- coldata$sample
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ cell_line + condition)

# Drop genes with almost no reads. They carry no information and slow things down.
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message("Genes kept after filtering: ", sum(keep))

dds <- DESeq(dds)
saveRDS(dds, snakemake@output[["dds"]])

# Shrink fold changes so genes with few reads do not get huge noisy values
res <- lfcShrink(dds, coef = "condition_dex_vs_control", type = "apeglm")

res_table <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by = "gene_id") %>%
  select(gene_id, gene_name, everything()) %>%
  arrange(padj)

write_tsv(res_table, snakemake@output[["results"]])

significant <- res_table %>%
  filter(!is.na(padj), padj < snakemake@params[["padj"]],
         abs(log2FoldChange) > snakemake@params[["lfc"]])
write_tsv(significant, snakemake@output[["significant"]])
message("Significant genes: ", nrow(significant),
        " (up: ", sum(significant$log2FoldChange > 0),
        ", down: ", sum(significant$log2FoldChange < 0), ")")

norm <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  rownames_to_column("gene_id") %>%
  left_join(gene_info, by = "gene_id") %>%
  select(gene_id, gene_name, everything())
write_tsv(norm, snakemake@output[["counts"]])
