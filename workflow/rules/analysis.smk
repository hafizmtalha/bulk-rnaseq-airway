rule deseq2:
    input:
        quant=expand("results/salmon/{sample}/quant.sf", sample=SAMPLES.index),
        tx2gene="data/reference/tx2gene.tsv",
        samples=config["samples"],
    output:
        dds="results/deseq2/dds.rds",
        results="results/deseq2/deseq2_results.tsv",
        significant="results/deseq2/significant_genes.tsv",
        counts="results/deseq2/normalized_counts.tsv",
    params:
        padj=config["padj_cutoff"],
        lfc=config["log2fc_cutoff"],
    log:
        "logs/deseq2.log",
    script:
        "../scripts/run_deseq2.R"


rule plots:
    input:
        dds="results/deseq2/dds.rds",
        results="results/deseq2/deseq2_results.tsv",
    output:
        pca="results/plots/pca.png",
        distances="results/plots/sample_distances.png",
        ma="results/plots/ma_plot.png",
        volcano="results/plots/volcano.png",
        heatmap="results/plots/top_genes_heatmap.png",
        genes="results/plots/top_genes_counts.png",
    params:
        padj=config["padj_cutoff"],
        lfc=config["log2fc_cutoff"],
    log:
        "logs/plots.log",
    script:
        "../scripts/make_plots.R"


rule gsea:
    input:
        results="results/deseq2/deseq2_results.tsv",
    output:
        table="results/gsea/hallmark_gsea.tsv",
        plot="results/plots/hallmark_gsea.png",
    log:
        "logs/gsea.log",
    script:
        "../scripts/run_gsea.R"
