# Gene set enrichment with the MSigDB hallmark sets.
# Genes are ranked by the DESeq2 Wald statistic, so we use every tested gene
# and not only the ones that passed a cutoff.

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)

res <- read_tsv(snakemake@input[["results"]], show_col_types = FALSE)

# lfcShrink drops the Wald statistic, so rebuild a ranking score from
# the sign of the fold change and the p value.
ranks <- res %>%
  filter(!is.na(pvalue), !is.na(gene_name)) %>%
  mutate(score = sign(log2FoldChange) * -log10(pvalue)) %>%
  group_by(gene_name) %>%
  slice_max(abs(score), n = 1, with_ties = FALSE) %>%
  ungroup()
ranks <- setNames(ranks$score, ranks$gene_name)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

gsea <- fgsea(pathways, ranks, minSize = 15, maxSize = 500)
gsea <- gsea %>%
  arrange(padj) %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ","))
write_tsv(gsea, snakemake@output[["table"]])

# Plot the 15 pathways with the smallest adjusted p values
plot_data <- gsea %>%
  slice_min(padj, n = 15, with_ties = FALSE) %>%
  mutate(pathway = gsub("_", " ", sub("HALLMARK_", "", pathway)),
         pathway = reorder(pathway, NES),
         direction = ifelse(NES > 0, "Up in dex", "Down in dex"),
         label = ifelse(padj < 0.05, "*", ""))

p <- ggplot(plot_data, aes(NES, pathway, fill = direction)) +
  geom_col() +
  geom_text(aes(label = label, hjust = ifelse(NES > 0, -0.5, 1.5)), size = 6, vjust = 0.75) +
  scale_fill_manual(values = c("Up in dex" = "#C44E52", "Down in dex" = "#4C72B0")) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "Normalized enrichment score", y = NULL, fill = NULL,
       title = "Hallmark pathways changed by dexamethasone",
       subtitle = "Top 15 pathways. * marks adjusted p value below 0.05") +
  theme_bw(base_size = 12)
ggsave(snakemake@output[["plot"]], p, width = 9, height = 6, dpi = 150)
