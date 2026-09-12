script 1:
uncov=$(bedtools coverage -sorted -a data/bedgraph/pilot.bed -b results/coverage/bedgraphs/22DX03-0786-201_bedgraph_nozero.bed -hist | grep all | head -n 1 | cut -f5)

script 2:
#this one OVERCOUNTS
total_bases=$(awk '{sum += $3-$2} END{print sum}' data/bedgraph/pilot.bed)

covered_bases=$(
    bedtools intersect \
        -sorted \
        -a data/bedgraph/pilot.bed \
        -b results/coverage/bedgraphs/22DX03-0786-201_bedgraph_nozero.bed \
        -wo |
    awk '{covered += $NF} END{print covered}'
)

uncov=$(awk -v c="$covered_bases" -v t="$total_bases" \
    'BEGIN{print 1-(c/t)}')

printf "22DX03-0786-201\t${uncov}\n" > results/coverage/uncov/samples/22DX03-0786-201_uncov.tsv

#merge the bedgraphs before estimate covered_bases
total_bases=$(awk '{sum += $3-$2} END{print sum}' data/bedgraph/pilot.bed)

#this return nothing, then uncov = 1
covered_bases=$(
    bedtools intersect \
        -a data/bedgraph/pilot.bed \
        -b <(bedtools sort -i results/coverage/bedgraphs/22DX03-0786-201_bedgraph_nozero.bed | bedtools merge -i -) \
        -wo |
    awk '{covered += $NF} END{print covered}'
)

uncov=$(awk -v c="$covered_bases" -v t="$total_bases" \
    'BEGIN{print 1-(c/t)}')

printf "22DX03-0786-201\t${uncov}\n" > results/coverage/uncov/samples/22DX03-0786-201_uncov.tsv

#use bedtools subtract - still the same 0.040059
# 1. Calculate total bases in your pilot.bed
total_bases=$(awk '{sum += $3-$2} END{print sum}' data/bedgraph/pilot.bed)

# 2. Directly calculate UNCOVERED bases using bedtools subtract
uncov_bases=$(
    bedtools subtract \
        -sorted \
        -a data/bedgraph/pilot.bed \
        -b results/coverage/bedgraphs/22DX03-0786-201_bedgraph_nozero.bed |
    awk '{uncov += $3-$2} END{print uncov}'
)

# 3. Calculate the fraction (uncov_bases / total_bases)
uncov=$(awk -v u="$uncov_bases" -v t="$total_bases" \
    'BEGIN{print u/t}')