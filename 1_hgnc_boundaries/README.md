# Get GENCODE-defined regions

> Specifically decided to use exons + slop rather than full transcripts as some HLA transcripts on ALT contigs at 10kb+
> in length. Moreover, LILAC specifically uses exons to select reads. The full HLA transcripts on the main contigs are
> fully covered by this approach, hence should capture all reads mapping elsewhere. However, we also explicitly add the
> HLA-ABC main contig transcripts with a 1000 bp slop to align with LILAC read slicing process.

Download required files
```bash
mkdir -p data/

urls='
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt
http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
'
parallel wget -P data/ ::: ${urls}
gzip -d data/gencode.v41.chr_patch_hapl_scaff.annotation.gff3.gz
```

Select all HLA genes, transcripts, and CDS across all hg38 contigs then merge regions
```bash
mkdir -p 1_regions/

for symbol in hla_{a,b,c}; do
  hgnc_symbol=$(tr '_' '-' <<< ${symbol^^});
  resp=$(
    curl \
        -s \
        -H 'Accept:application/json' \
        "https://rest.genenames.org/search/symbol/${hgnc_symbol}"
  );
  hgnc_id=$(jq -r <<< ${resp} '.response.docs[] | .hgnc_id');
  echo ${symbol} ${hgnc_id};
done | tee 1_regions/hgnc_ids.txt

./scripts/get_regions.py > 1_regions/regions_gencode_unmerged.bed
./scripts/convert_gencode_ucsc_contigs.py | sort | uniq > 1_regions/regions_ucsc_unmerged.bed

bedtools sort -i 1_regions/regions_ucsc_unmerged.bed | \
  bedtools slop -b 500 -g data/hg38.chrom.sizes -i - | \
  bedtools merge -c 4 -o distinct -i - > 1_regions/regions_ucsc.bed

# Add in full main HLA-ABC genes with 1000 bp slop to mirror HMF read selection
grep '^chr6.\+\tgene\t.\+gene_name=HLA-[ABC]' data/gencode.v41.chr_patch_hapl_scaff.annotation.gff3 | \
  awk '{ print $1, $4-1, $5 }' | \
  tr ' ' '\t' | \
  bedtools sort -i - | \
  bedtools slop -b 1000 -g data/hg38.chrom.sizes -i - > 1_regions/regions_main_genes.bed

# Combine
cat \
  1_regions/regions_main_genes.bed \
  <(cut -f1-3 -d$'\t' 1_regions/regions_ucsc.bed) | \
  bedtools sort -i - | \
  bedtools merge -i - > regions.bed
```
