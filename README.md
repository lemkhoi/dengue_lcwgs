#Dengue low-coverage WGS workflow

This directory is a working snapshot of a Snakemake workflow for low-coverage whole-genome sequencing (LC-WGS) of HCV and dengue cohorts. It trims paired FASTQ reads, aligns them to GRCh38, removes duplicates, calculates quality-control (QC) metrics, performs QUILT genotype imputation, filters imputed variants, and optionally performs HLA calling/imputation and downstream analyses.

The workflow was adapted from [LCWGS_pipeline](https://github.com/Suuuuuuuus/LCWGS_pipeline). The local version is project-specific and contains absolute BMRC paths, exploratory scripts, and some modules that are not connected to the current master workflow.

## Main workflow

Run targets from `pipelines/master_lc.smk` in this order after a successful dry-run:

1. **Alignment and QC:** `alignment_and_qc_all` trims adapters, aligns to GRCh38, fixes read metadata, removes duplicates, indexes BAMs, calculates coverage/duplication/sex-chromosome metrics, and writes `results/metrics.tsv` plus QC plots under `graphs/qc/`.
2. **Genome-wide imputation:** `imp_all` prepares regional QUILT references, imputes each batch, concatenates regional VCFs, merges batches, and applies allele-frequency, INFO-score, and Omni5M-site filters. Principal outputs are under `results/imputation/` and `results/wip_vcfs/`.
3. **HLA analysis:** `hla_all` runs HLA-LA, builds HLA alignment matrices, prepares the HLA reference, and performs QUILT-HLA imputation. Outputs are under `results/hla/`.
4. **Optional targets:** `qc_all` reruns QC without requiring alignment; `sv_all` and `hla_coverage_all` generate coverage inputs for structural-variant/HLA-region work.
5. **Second-stage TOPMed imputation:** use only carefully reviewed VCFs from `results/wip_vcfs/`. The scripts under `scripts/` show a previous manual submission and post-processing route, not a portable production workflow.

## Important cautions

- This is **not a self-contained release**: `data/`, scheduler wrappers, environment files, and several scripts referenced by standalone modules are absent.
- Paths and sample IDs are a mixture of the original handover and later edits. Search for `/well/ansari/` before running: `grep -R "/well/ansari/" pipelines scripts`.
- `scripts/1.8_submit_topmed.sh` contains an embedded API token. Treat it as compromised: revoke/rotate it, remove it from version history, and load a replacement from a protected environment variable. Do not run or commit the file as written.
- Several files in `scripts/` are one-off diagnostics. In particular, `1.3_uncov_22DX03-0786-201.sh` is an experiment rather than valid end-to-end shell code, and `1.4_sex_check_combined.sh` currently requests three fields from `bcftools query` while its AWK code expects four (including allele depth).
- Check genome-build and chromosome naming consistency (`38`/`hg38`, `chr1` versus `1`) at every hand-off.


### `pipelines/`

| File | Description |
|---|---|
| `alignment.smk` | Aligns cleaned paired reads with BWA, sorts/indexes BAMs, adds read groups, fixes mates, removes duplicates with Picard, and builds batch BAM lists. |
| `auxiliary.smk` | Shared Python helpers for reading sample TSVs and expanding chromosome-region VCF, BAM, and coverage paths. |
| `config.json` | Project configuration: paths, batches, sample lists, reference panels, QUILT settings, QC thresholds, HLA coordinates/version, and GWAS settings. |
| `coverage.smk` | Calculates genome-wide, HLA-region, uncovered-fraction, and sex-chromosome coverage metrics and aggregates them across samples. |
| `dup_rate.smk` | Extracts duplicate fractions from Picard metrics and combines per-sample values. |
| `fastqc.smk` | Runs FastQC/MultiQC and derives sequencing-error/BQSR-related summaries; it is not included by the current master workflow. |
| `filter_vcf.smk` | Merges batch-level QUILT VCFs and filters variants by gnomAD MAF, QCTool INFO score, Omni5M sites, and HLA scaffold sites. |
| `hla.smk` | Runs HLA-LA on each BAM and writes per-sample HLA calls. |
| `hla_imputation_method.smk` | Prepares a QUILT-HLA reference, calculates allele-alignment matrices, and imputes HLA-A/B/C/DRB1/DQB1. |
| `hla_ref_panel.smk` | Builds an HLA-region reference by merging GAMCC and 1000 Genomes data in chunks; standalone/not included by `master_lc.smk`. |
| `imputation.smk` | Prepares regional QUILT reference objects, imputes each batch/region, and concatenates regional VCFs by chromosome. |
| `imputation_prep.smk` | Prepares BAM lists, recombination maps, reference haplotype/legend files, and imputation chunks; currently commented out in the master. |
| `master_lc.smk` | Main entry point defining `alignment_and_qc_all`, `qc_all`, `imp_all`, `hla_all`, `sv_all`, and `hla_coverage_all`. |
| `merge_boost.smk` | One-off rules that merge FASTQs from repeated dengue sequencing batches. |
| `nonahore.smk` | Simulation and structural-variant evaluation targets for Nonahore; references Python scripts not included in this snapshot. |
| `phasing.smk` | Prepares/phases GAMCC and 1000 Genomes HLA alleles with Beagle and assesses phasing concordance. |
| `population.smk` | Converts chip and LC-WGS data to BGEN and performs PCA/projection analyses for population structure. |
| `post_gw.smk` | Lifts imputation-server genome-wide VCF output between genome builds. |
| `post_hla.smk` | Lifts HLA scaffold data and creates chromosome 6 inputs for two- or three-stage HLA imputation. |
| `preprocess.smk` | Trims adapters from paired FASTQs with Trimmomatic and writes temporary cleaned read pairs. |
| `qc.smk` | Merges sample metrics and creates overall/per-batch coverage, HLA coverage, QC, and sex-coverage plots. |
| `reference.smk` | Builds/indexes GRCh38 and constructs chunked merged reference panels from 1000 Genomes and MalariaGEN data. |
| `regions.smk` | Keeps primary alignments and calculates chunked whole-genome/HLA coverage with Coverotron. |
| `software.smk` | Central mapping of external executable/JAR paths, including QCTool, Picard, QUILT-HLA, HLA-LA, IMPUTE2, Beagle, and coverage tools. |
| `test.smk` | Small test workflow for generating a nonzero-coverage bedGraph for one BAM. |
| `vanilla_alignment.yml` | Minimal Conda environment for Snakemake, BWA, samtools, Picard/Picard Slim, and Java 17. |

### `scripts/`

| File | Description |
|---|---|
| `1.0_metrics_check.ipynb` | Exploratory notebook that merges QC metrics, identifies missing results, examines coverage/duplication distributions, and investigates sex discordance. |
| `1.0_plot_functions.ipynb` | Exploratory notebook for QUILT VCF inspection and INFO-score versus allele-frequency summaries/plots. |
| `1.1_test_vanilla_aln.smk` | Single-sample test of the BWA/samtools/Picard alignment and duplicate-removal procedure. |
| `1.1_test_vanilla_aln_all.sh` | Slurm script that applies the same alignment procedure to a hard-coded six-sample test set. |
| `1.2_picard.sh` | One-sample Slurm diagnostic that reruns Picard duplicate marking on an existing BAM. |
| `1.3_dup_metric_combined.sh` | AWK utility that combines Picard `DuplicationMetrics` files into one TSV. |
| `1.3_uncov_22DX03-0786-201.sh` | Scratch comparison of three Bedtools approaches for one sample's uncovered fraction; contains literal section labels and should not be run unchanged. |
| `1.4_sex_check.sh` | Slurm-array workflow that calls non-PAR chromosome X variants per BAM and runs PLINK `--check-sex`. |
| `1.4_sex_check_combined.sh` | Aggregates PLINK sex-check output with heterozygosity, depth, allele-balance, and genotype-quality summaries; currently needs its `bcftools query` field list fixed. |
| `1.4_sex_guess_ploidy.sh` | Calls chromosome X non-PAR variants and runs the bcftools `+guess-ploidy` plugin for listed samples. |
| `1.5_vis_imputation_plots4.smk` | Pandas/`lcwgsus` workflow that merges imputed-variant and allele-frequency tables and plots INFO score by allele-frequency bin. |
| `1.5_vis_imputation_plots_bcf.smk` | Bcftools-based alternative that annotates allele frequencies, extracts INFO fields, combines chromosomes, and plots INFO score by frequency bin. |
| `1.6_filter_vcf_manual.sh` | Manual QCTool route converting chromosome VCFs to BGEN and calculating SNP INFO statistics. |
| `1.7_filter_samples_cov40.sh` | Keeps a supplied sample list, filters variants to `INFO_SCORE > 0.9`, and indexes chromosome VCFs for TOPMed submission. |
| `1.8_submit_topmed.sh` | Submits filtered VCFs to the TOPMed imputation API; unsafe as stored because it embeds a credential. |
| `1.9_prepare_plink.smk` | Filters returned TOPMed VCFs to `R2 > 0.5` and `MAF > 0.05`, then concatenates autosomes; despite its name, it does not create PLINK files. |

## Expected directory layout

```text
data/
  fastq/                 # paired input FASTQs
  references/            # GRCh38 FASTA plus indexes/dictionary
  sample_tsvs/           # one-column sample lists and per-batch lists
  bedgraph/               # accessible-genome and HLA BED files
  imputation_accessories/# chunk definitions and recombination maps
  ref_panel/              # imputation reference data
results/                  # BAMs, metrics, imputation, filtered VCFs, HLA outputs
graphs/                   # QC and imputation plots
pipelines/                # Snakemake workflow modules
scripts/                  # diagnostics and one-off downstream analyses
```

Sample TSVs are expected to be one column with no header unless a specific script says otherwise. Create output directories through Snakemake where possible rather than manually moving partial results.
