rule fastp:
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_1.fastq.gz",
        r2="results/trimmed/{sample}_2.fastq.gz",
        html="results/fastp/{sample}.html",
        json="results/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    threads: config["threads"]
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1} --out2 {output.r2} \
              --detect_adapter_for_pe --thread {threads} \
              --html {output.html} --json {output.json} 2> {log}
        """


rule multiqc:
    input:
        expand("results/fastp/{sample}.json", sample=SAMPLES.index),
        expand("results/salmon/{sample}/quant.sf", sample=SAMPLES.index),
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        "multiqc --force --outdir results/multiqc results/fastp results/salmon 2> {log}"
