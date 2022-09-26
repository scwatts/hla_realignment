# HLA realignment

## Background

LILAC requires that reads relevant to HLA-C, HLA-B, and HLA-C are aligned to these genes on the main chr6 contig and
will not consider read alignments on ALT contigs. Hence, to use LILAC with the hg38 assembly all reads mapped to ALT
sequences related to these HLA genes must be identified and realigned to the main chr6 contig.

Here, I present an approach that uses a set of definied regions to collect reads relevant to HLA-A, HLA-B, and HLA-C on
ALT contigs and realign these reads appropriately. The major steps are:

* define a set of regions
  * GENCODE exons
  * homologoues sequences
* use the defined set of regions to extract reads
* align extracted reads to hg38 chr6
* run LILAC


## Caveats

* alignment not done directly to HLA-A, HLA-B, HLA-C since LILAC requires coordinates in context of chr6
* alignment restricted to chr6 with the aim to save time (not assessed); better alignments may exist in other contigs
* slicing can generate orphan reads and these are aligned as SE, undetermined how LILAC processes SE reads
* LILAC can use PURPLE CNV data and SAGE variant calls but these would be generated using different alignments
