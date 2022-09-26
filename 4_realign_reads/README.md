## Collect and map reads

> We do not search the BAM for read mates not within target regions since LILAC filters these anyway, see:
> https://github.com/hartwigmedical/hmftools/blob/b49d428/lilac/src/main/java/com/hartwig/hmftools/lilac/read/BamRecordReader.java#L300-L325


```bash
gds_base_dir=gds://production/analysis_data/SBJ02737/wgs_tumor_normal/2022082156da5dce

bam_gds_tumor=${gds_base_dir}/L2201206_L2201205_dragen/PRJ221999_tumor.bam
bam_gds_normal=${gds_base_dir}/L2201206_L2201205_dragen/PRJ221998_normal.bam

# Get chr6 from hg38 so that alignment coords are compatible with LILAC
mkdir -p data/
samtools faidx \
  -o data/ref_extracted.fa \
  ~/repos/reference_data/reference_data/genomes/hg38.fa \
  chr6
bwa index data/ref_extracted.fa

# Get reads for HLA tumor and normal and realign to hg38 chr6
./scripts/process_reads.sh ${bam_gds_tumor} data/ref_extracted.fa
./scripts/process_reads.sh ${bam_gds_normal} data/ref_extracted.fa
```
