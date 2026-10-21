"""Build a transcript to gene table from the GENCODE FASTA headers.

A GENCODE header looks like this:
>ENST00000456328.2|ENSG00000290825.1|...|DDX11L2-202|DDX11L2|1657|lncRNA|

Salmon run with --gencode keeps only the first field as the transcript name,
so we keep that and pair it with the gene id, the gene name and the
transcript type.
"""

import gzip

fasta = snakemake.input[0]
out = snakemake.output[0]

with gzip.open(fasta, "rt") as f, open(out, "w") as o:
    o.write("transcript_id\tgene_id\tgene_name\ttranscript_type\n")
    for line in f:
        if not line.startswith(">"):
            continue
        fields = line[1:].strip().split("|")
        o.write(f"{fields[0]}\t{fields[1]}\t{fields[5]}\t{fields[7]}\n")
