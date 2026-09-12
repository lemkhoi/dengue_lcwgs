configfile: "pipelines/config.json"
home_dir = config['home_dir']

import pandas as pd
samples_lc = read_tsv_as_lst(config['samples_lc'])

rule extract_picard_dup_rate:
    input:
        file = 'data/bams/{id}.metrics.txt'
    output:
        duprate = "results/dup_rate/samples/dup_rate_{id}.tsv"
    localrule: True
    shell: """
        mkdir -p results/dup_rate/samples/

        duprate=$(sed -n '8p' {input.file} | cut -f9)
        echo "{wildcards.id}\t$duprate" >> {output.duprate}
    """

rule aggregate_dup_rate:
    input:
        files = expand("results/dup_rate/samples/dup_rate_{id}.tsv", id = samples_lc)
    output:
        duprate = "results/dup_rate/dup_rate_all.tsv"
    shell: """
        cat {input.files} >> {output.duprate}
    """

