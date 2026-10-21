rule salmon_index:
    input:
        "data/reference/transcripts.fa.gz",
    output:
        directory("data/reference/salmon_index"),
    log:
        "logs/salmon_index.log",
    threads: config["threads"]
    shell:
        "salmon index --gencode -t {input} -i {output} -p {threads} > {log} 2>&1"


rule tx2gene:
    input:
        "data/reference/transcripts.fa.gz",
    output:
        "data/reference/tx2gene.tsv",
    script:
        "../scripts/make_tx2gene.py"


rule salmon_quant:
    input:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz",
        index="data/reference/salmon_index",
    output:
        "results/salmon/{sample}/quant.sf",
    log:
        "logs/salmon/{sample}.log",
    threads: config["threads"]
    shell:
        """
        salmon quant -i {input.index} -l A \
            -1 {input.r1} -2 {input.r2} \
            --validateMappings --gcBias --seqBias \
            -p {threads} -o results/salmon/{wildcards.sample} > {log} 2>&1
        """
