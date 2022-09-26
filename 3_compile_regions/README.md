# Merge regions

Require regions to be at least 50 bp. Add in all HLA contigs.

```bash
hla_contigs=$(
  awk < ~/repos/reference_data/reference_data/genomes/hg38.fa.fai \
    '
      BEGIN { OFS="\t" }
      /^HLA/ { print $1, "0", $2 }
    '
)

cat \
 ../1_hgnc_boundaries/regions.bed \
 <(cut -f1-3 -d$'\t' ../2_homologous_regions/regions.bed) \
 <(echo "${hla_contigs}") | \
 bedtools sort -i - | \
 bedtools merge -i - | \
 awk '{ s=($3+1)- $2; if (s >= 50) { print } }' > regions.bed
```
