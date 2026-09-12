configfile: "pipelines/config.json"
include: "software.smk"
include: "auxiliary.smk"
home_dir = config['home_dir']

rule vanilla_alignment:
    input:
        fastq1 = "data/fastq_cleaned/{id}_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/{id}_2.fastq.gz",
        reference = "data/references/GRCh38.fa"
    output:
        bam = "data/bams/{id}.bam",
        bai = "data/bams/{id}.bam.bai",
        tmp1 = temp("data/bams/{id}.tmp1.bam"),
        metric = "data/bams/{id}.metrics.txt"
    resources:
        mem = '40G'
    params: 
        sample = "{id}",
        picard = tools["picard_plus"]
    threads: 6
    shell: """
        bwa mem -t {threads} {input.reference} {input.fastq1} {input.fastq2} | samtools view -b -o {output.tmp1}
        
        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}

        samtools index {output.bam}

        picard AddOrReplaceReadGroups \
        -VERBOSITY ERROR \
        -I {output.bam} \
        -O {output.tmp1} \
        -RGLB OGC \
        -RGPL Illumina \
        -RGPU unknown \
        -RGSM {params.sample}

        {params.picard} FixMateInformation -I {output.tmp1}

        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}

        {params.picard} MarkDuplicates \
        -I {output.bam} \
        -O {output.tmp1} \
        -M {output.metric} \
        --REMOVE_DUPLICATES

        samtools sort -@6 -m 1G -o {output.bam} {output.tmp1}
        samtools index {output.bam}
    """

def get_bams_from_batch(wildcards):
    file = f"data/sample_tsvs/batches/samples_{wildcards.batch}.tsv"
    return [f"data/bams/{s}.bam" for s in read_tsv_as_lst(file)]

rule prepare_bamlist:
    input:
        bams = get_bams_from_batch
    output:
        bamlist = "results/imputation/bamlists/{batch}/bamlist.txt"
    localrule: True
    shell: """
        mkdir -p results/imputation/bamlists/{wildcards.batch}/

        ls {input.bams} > {output.bamlist}
    """