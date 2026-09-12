home_dir= "/well/ansari/users/osr869/"
rule all:
    input:
        bam = "data/bams/28EI-003-084.bam",
        bai = "data/bams/28EI-003-084.bam.bai",
        metric = "data/bams/28EI-003-084.metrics.txt"

rule vanilla_alignment:
    input:
        fastq1 = "data/fastq_cleaned/28EI-003-084_1.fastq.gz",
        fastq2 = "data/fastq_cleaned/28EI-003-084_2.fastq.gz",
        reference = "data/references/GRCh38.fa"
    output:
        bam = "data/bams/28EI-003-084.bam",
        bai = "data/bams/28EI-003-084.bam.bai",
        tmp1 = temp("data/bams/28EI-003-084.tmp1.bam"),
        metric = "data/bams/28EI-003-084.metrics.txt"
    resources:
        mem = '40G'
    params: 
        sample = "28EI-003-084",
        picard = f"java -Xmx40G -Xms20G -jar {home_dir}conda/skylake/envs/vani_aln17/share/picard-slim-2.27.4-0/picard.jar"
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