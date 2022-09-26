# Get homologous regions

Done under the assumption that we don't really care about reads mapped to non-homologous HLA regions since these should
never realign back to the HLA genes on the main contigs.

Search subject

```bash
mkdir -p 1_homologies/

contigs=$(awk '$1 ~ /chr6|chrUn/ { print $1 }' ~/repos/reference_data/reference_data/genomes/hg38.fa.fai)

samtools faidx > 1_homologies/genome_contigs.fasta \
  ~/repos/reference_data/reference_data/genomes/hg38.fa ${contigs}
```

Search query, rename FASTA ids

```bash
awk '{ print $1 ":" $2 "-" $3 }' ~/repos/reference_data/reference_data/hmftools/lilac/hla.38.bed |
    samtools faidx > 1_homologies/hmf_hla_genes.fasta \
    ~/repos/reference_data/reference_data/genomes/hg38.fa \
    --region-file -

rename_map=$(cat <<EOF | tr ' ' '\t'
chr6:29940260-29946884 HLA-A
chr6:31352872-31358188 HLA-B
chr6:31267749-31273130 HLA-C
EOF
)
tmp_fp=$(mktemp $(pwd -P)/tmp.XXXXXX)
seqkit replace > ${tmp_fp} \
  -p '(.+)$' \
  -r '{kv}' \
  -k <(echo "${rename_map}") \
  1_homologies/hmf_hla_genes.fasta
mv ${tmp_fp} 1_homologies/hmf_hla_genes.fasta
```

Run search

```bash
blat -out=blast9 1_homologies/genome_contigs.fasta 1_homologies/hmf_hla_genes.fasta 1_homologies/results.tsv
```

Processing results

```bash
grep -v '#' 1_homologies/results.tsv | \
  awk '
    BEGIN {
      OFS="\t";
    }
    $4 >= 50 {
      d[1] = $9;
      d[2] = $10;
      asort(d);
      print $2, d[1]-1, d[2], "gene:" $1 ";type=blat_homologous"
    }
  ' | \
  bedtools sort -i - | \
  bedtools merge -c 4 -o distinct -i - > regions.bed
```
