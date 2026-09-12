awk '
BEGIN { OFS="\t" }

FNR == 1 {
    sample = FILENAME
    sub(/^.*\//, "", sample)
    sub(/\.metrics\.txt$/, "", sample)
}

$0 ~ /^## METRICS CLASS[[:space:]]+picard\.sam\.DuplicationMetrics/ {
    getline header
    getline data

    if (!printed_header) {
        print "SAMPLE", header
        printed_header = 1
    }

    print sample, data
}
' data/bams/*.metrics.txt > ../results/duplication_metrics_combined.tsv