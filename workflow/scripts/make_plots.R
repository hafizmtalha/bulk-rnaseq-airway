# Make the main figures: sample level checks first, then gene level results.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
  library(tidyr)
})

theme_set(theme_bw(base_size = 13))
padj_cut <- snakemake@params[["padj"]]
lfc_cut <- snakemake@params[["lfc"]]

dds <- readRDS(snakemake@input[["dds"]])
res <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

# The variance stabilising transform puts counts on a log like scale
# where every gene has similar spread. Good for PCA and clustering.
vsd <- vst(dds, blind = FALSE)

# PCA of samples
pca <- plotPCA(vsd, intgroup = c("condition", "cell_line"), returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, shape = cell_line)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1 (", pct[1], "% of variance)"),
       y = paste0("PC2 (", pct[2], "% of variance)"),
       title = "Samples split by treatment on PC1",
       colour = "Condition", shape = "Cell line") +
  scale_colour_manual(values = c(control = "#4C72B0", dex = "#DD8452"))
ggsave(snakemake@output[["pca"]], p, width = 7, height = 5, dpi = 150)

# Sample to sample distances
dists <- dist(t(assay(vsd)))
mat <- as.matrix(dists)
labels <- paste(vsd$cell_line, vsd$condition, sep = " ")
rownames(mat) <- labels
colnames(mat) <- labels
png(snakemake@output[["distances"]], width = 1000, height = 850, res = 150)
pheatmap(mat, clustering_distance_rows = dists, clustering_distance_cols = dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Distance between samples")
dev.off()

# Label each gene as up, down or not changed
res <- res %>%
  filter(!is.na(padj)) %>%
  mutate(status = case_when(
    padj < padj_cut & log2FoldChange > lfc_cut ~ "Up in dex",
    padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down in dex",
    TRUE ~ "Not changed"
  ))
status_cols <- c("Up in dex" = "#C44E52", "Down in dex" = "#4C72B0", "Not changed" = "grey75")
top <- res %>% filter(status != "Not changed") %>% slice_min(padj, n = 15)

# MA plot
p <- ggplot(res, aes(log10(baseMean), log2FoldChange, colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_colour_manual(values = status_cols) +
  labs(x = "Mean expression (log10 of normalized counts)",
       y = "log2 fold change (dex vs control)",
       title = "MA plot", colour = NULL)
ggsave(snakemake@output[["ma"]], p, width = 8, height = 5, dpi = 150)

# Volcano plot
p <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = top, aes(label = gene_name), colour = "black",
                  size = 3.5, max.overlaps = 30) +
  scale_colour_manual(values = status_cols) +
  labs(x = "log2 fold change (dex vs control)",
       y = "-log10 adjusted p value",
       title = paste0(sum(res$status == "Up in dex"), " genes up and ",
                      sum(res$status == "Down in dex"), " genes down after dexamethasone"),
       colour = NULL)
ggsave(snakemake@output[["volcano"]], p, width = 8, height = 6, dpi = 150)

# Heatmap of the 40 most significant genes
top40 <- res %>% filter(status != "Not changed") %>% slice_min(padj, n = 40)
mat <- assay(vsd)[top40$gene_id, ]
rownames(mat) <- top40$gene_name
mat <- mat - rowMeans(mat)
annot <- data.frame(Condition = vsd$condition, `Cell line` = vsd$cell_line,
                    row.names = colnames(vsd), check.names = FALSE)
png(snakemake@output[["heatmap"]], width = 1000, height = 1400, res = 150)
pheatmap(mat, annotation_col = annot, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = "Top 40 genes (centered expression)")
dev.off()

# Normalized counts for the six most significant genes, with lines joining
# the treated and control sample from the same cell line
norm <- counts(dds, normalized = TRUE)
top6 <- res %>% slice_min(padj, n = 6)
long <- as.data.frame(t(norm[top6$gene_id, ]))
colnames(long) <- top6$gene_name
long$condition <- dds$condition
long$cell_line <- dds$cell_line
long <- pivot_longer(long, cols = all_of(top6$gene_name), names_to = "gene", values_to = "count")
long$gene <- factor(long$gene, levels = top6$gene_name)
p <- ggplot(long, aes(condition, count)) +
  geom_line(aes(group = cell_line), colour = "grey60") +
  geom_point(aes(colour = condition), size = 3) +
  scale_y_log10() +
  facet_wrap(~gene, scales = "free_y") +
  scale_colour_manual(values = c(control = "#4C72B0", dex = "#DD8452")) +
  labs(x = NULL, y = "Normalized count (log scale)",
       title = "Most significant genes, one line per cell line") +
  theme(legend.position = "none")
ggsave(snakemake@output[["genes"]], p, width = 9, height = 6, dpi = 150)
