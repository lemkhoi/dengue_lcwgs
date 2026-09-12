#!/bin/bash

TOKEN="eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJsZW1raG9pIiwibmJmIjoxNzgxMjc4NzYxLCJtYWlsIjoibGluYTQ3NzJAb3guYWMudWsiLCJhcGlfaGFzaCI6IlFQbUlvVlpWUGpnWnNMV2NWT0hvRWdOc1dYcGNyciIsInJvbGVzIjpbXSwiaXNzIjoiY2xvdWRnZW5lIiwibmFtZSI6Iktob2kgTGUiLCJhcGkiOnRydWUsImV4cCI6MTc4Mzg3MDc2MSwidG9rZW5fdHlwZSI6IkFQSV9UT0tFTiIsImlhdCI6MTc4MTI3ODc2MSwidXNlcm5hbWUiOiJsZW1raG9pIn0.Cr2I5hirN7yF0Oy40SIJwxEAyF-v1YOl87P3cUdDT9A"

filtered_dir="/well/ansari/users/osr869/host/lcwgs_edits/topmed/covered40_samples"
FILES=(${filtered_dir}/*.vcf.gz);
ARGS=();
for f in "${FILES[@]}"; do
    ARGS+=(-F "files=@$f")
done;
curl -H "X-Auth-Token: ${TOKEN}" \
    "${ARGS[@]}" \
    -F "refpanel=topmed-r3" \
    -F "build=hg38" \
    -F "phasing=eagle" \
    -F "population=all" \
    -F "meta=yes" \
    'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2/jobs/submit/imputationserver2'

