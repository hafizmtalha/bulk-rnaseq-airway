def fastq_url(wildcards):
    row = SAMPLES.loc[wildcards.sample]
    return row["fq1"] if wildcards.read == "1" else row["fq2"]


rule download_fastq:
    output:
        "data/fastq/{sample}_{read}.fastq.gz",
    params:
        url=fastq_url,
        test_mode=config["test_mode"],
        lines=config["test_read_pairs"] * 4,
    log:
        "logs/download/{sample}_{read}.log",
    wildcard_constraints:
        read="1|2",
    retries: 2
    shell:
        """
        if [ "{params.test_mode}" = "True" ]; then
            # Stream the file and stop after the first N reads.
            # curl complains when the stream is cut early, so we ignore that.
            set +o pipefail
            curl -sL {params.url} | gzip -dc 2> /dev/null | head -n {params.lines} | gzip > {output} 2> {log}
        else
            # Big files can drop halfway. Keep a partial file and resume it.
            for attempt in 1 2 3 4 5 6 7 8 9 10; do
                if curl -L --retry 5 --retry-all-errors -C - -o {output}.part {params.url} 2>> {log}; then
                    break
                fi
                sleep 30
            done
            mv {output}.part {output}
        fi
        gzip -t {output}
        """


rule download_transcriptome:
    output:
        "data/reference/transcripts.fa.gz",
    params:
        url=config["transcriptome_url"],
    log:
        "logs/download/transcriptome.log",
    shell:
        "curl -L --retry 5 -o {output} {params.url} 2> {log}"
