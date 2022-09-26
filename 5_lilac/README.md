```bash
docker run -ti \
  -v ~/repos/reference_data/reference_data/:/reference_data/ \
  -v $(pwd -P)/../:/working/ \
  -w /working/5_lilac \
  scwatts/lilac:1.2

java \
  -Xmx4g \
  -jar /opt/lilac/lilac.jar \
    -sample SBJ02737 \
    -reference_bam ../4_realign_reads/2_alignments/PRJ221998_normal/PRJ221998_normal_realigned.bam \
    -tumor_bam ../4_realign_reads/2_alignments/PRJ221999_tumor/PRJ221999_tumor_realigned.bam \
    -ref_genome_version 38 \
    -ref_genome /reference_data/genomes/hg38.fa \
    -resource_dir /reference_data/hmftools/lilac/ \
    -threads 1 \
    -write_all_files \
    -output_dir 1_output/
```
