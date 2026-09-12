#!/usr/bin/env bash
set -euo pipefail

outdir="/well/ansari/users/osr869/host/lcwgs_edits/results/X_check"
summary_out="${outdir}/sexcheck_summary.tsv"

printf "sample_id\tFID\tIID\tPEDSEX\tSNPSEX\tSTATUS\tF\tHET_COUNT\tCALLED_SITES\tHET_RATE\tMEAN_DP\tALLELE_BALANCE_MEAN\tGQ_0_9\tGQ_10_19\tGQ_20_29\tGQ_30_39\tGQ_40_PLUS\tGQ_MISSING\n" > "$summary_out"

for sexfile in "${outdir}"/*_sexcheck.sexcheck; do
    [[ -e "$sexfile" ]] || continue

    sample_id=$(basename "$sexfile" _sexcheck.sexcheck)
    bcf="${outdir}/${sample_id}_chrX_nonPAR.bcf"

    if [[ ! -f "$bcf" ]]; then
        echo "Missing BCF for ${sample_id}: ${bcf}" >&2
        continue
    fi

    stats=$(
        bcftools query -f '[%GT\t%DP\t%GQ\n]' "$bcf" |
        awk -F'\t' '
            function norm(gt) {
                gsub(/\|/, "/", gt)
                return gt
            }

            function is_missing_gt(gt,   n, a, i) {
                gt = norm(gt)
                if (gt == "." || gt == "./." || gt == ".") return 1
                n = split(gt, a, "/")
                for (i = 1; i <= n; i++) {
                    if (a[i] == ".") return 1
                }
                return 0
            }

            function is_het_gt(gt,   n, a, i) {
                gt = norm(gt)
                if (gt == "." || gt == "./." || gt == ".") return 0
                n = split(gt, a, "/")
                for (i = 2; i <= n; i++) {
                    if (a[i] != a[1]) return 1
                }
                return 0
            }

            BEGIN {
                het=0
                called=0
                sumdp=0
                ndp=0
                sumab=0
                nab=0

                gq0_9=0
                gq10_19=0
                gq20_29=0
                gq30_39=0
                gq40p=0
                gqmiss=0
            }

            {
                gt = $1
                dp = $2
                ad = $3
                gq = $4

                if (!is_missing_gt(gt)) {
                    called++

                    if (dp != "." && dp != "") {
                        sumdp += dp
                        ndp++
                    }

                    if (is_het_gt(gt)) {
                        het++

                        # Allele balance from AD at heterozygous sites:
                        # AB = alt_depth / (ref_depth + sum(alt_depths))
                        if (ad != "." && ad != "") {
                            n = split(ad, a, ",")
                            if (n >= 2) {
                                ref = a[1] + 0
                                alt = 0
                                for (i = 2; i <= n; i++) {
                                    if (a[i] != ".") alt += a[i] + 0
                                }
                                denom = ref + alt
                                if (denom > 0) {
                                    sumab += alt / denom
                                    nab++
                                }
                            }
                        }
                    }
                }

                if (gq == "." || gq == "") {
                    gqmiss++
                } else {
                    gq = gq + 0
                    if (gq < 10) gq0_9++
                    else if (gq < 20) gq10_19++
                    else if (gq < 30) gq20_29++
                    else if (gq < 40) gq30_39++
                    else gq40p++
                }
            }

            END {
                het_rate = (called > 0) ? het / called : 0
                mean_dp = (ndp > 0) ? sumdp / ndp : "NA"
                mean_ab = (nab > 0) ? sumab / nab : "NA"

                printf "%d\t%d\t%.6f\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    het, called, het_rate, mean_dp, mean_ab,
                    gq0_9, gq10_19, gq20_29, gq30_39, gq40p, gqmiss
            }
        '
    )

    IFS=$'\t' read -r het_count called_sites het_rate mean_dp mean_ab gq0_9 gq10_19 gq20_29 gq30_39 gq40p gqmiss <<< "$stats"

    awk -v sample="$sample_id" \
        -v het="$het_count" \
        -v called="$called_sites" \
        -v rate="$het_rate" \
        -v mdp="$mean_dp" \
        -v mab="$mean_ab" \
        -v gq0_9="$gq0_9" \
        -v gq10_19="$gq10_19" \
        -v gq20_29="$gq20_29" \
        -v gq30_39="$gq30_39" \
        -v gq40p="$gq40p" \
        -v gqmiss="$gqmiss" '
        NR > 1 {
            print sample "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" het "\t" called "\t" rate "\t" mdp "\t" mab "\t" gq0_9 "\t" gq10_19 "\t" gq20_29 "\t" gq30_39 "\t" gq40p "\t" gqmiss
        }
    ' "$sexfile" >> "$summary_out"
done

echo "$summary_out"