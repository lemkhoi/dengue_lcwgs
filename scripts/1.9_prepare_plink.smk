home_dir = "/well/ansari/users/osr869/host/lcwgs_edits/"
chrs = [i for i in range(1,23)]

rule all:
    input:
        home_dir + "gwas/wg.topmed.vcf.gz"

rule filter_topmed_vcf:
    input:
        vcfs = home_dir + "topmed/covered40_result/chr{chr}.dose.vcf.gz"
    output:
        vcf_filtered = home_dir + "topmed/covered40_result_filtered/topmed.filtered.chr{chr}.vcf.gz"
    shell: """
        bcftools view \
            -i 'INFO/R2 > 0.5 && INFO/MAF > 0.05' \
            -c 1 \
            {input.vcfs} -Oz -o {output.vcf_filtered}
        tabix -f {output.vcf_filtered}
    """

rule merge_wg_topmed_vcf:
    input:
        vcf_filtered = expand(home_dir + "topmed/covered40_result_filtered/topmed.filtered.chr{chr}.vcf.gz", chr=chrs)
    output:
        vcf_wg = home_dir + "gwas/wg.topmed.vcf.gz"
    shell: """
        bcftools concat -Oz -o {output.vcf_wg} {input.vcf_filtered}
        tabix -f {output.vcf_wg}
    """