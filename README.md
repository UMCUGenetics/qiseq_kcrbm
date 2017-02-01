# Introduction
R code for Kernel Convolved Rule Based Mapping (KC-RBM) as described by De Jong et al. KC-RBM exploits distance, orientation and insertion density across tumors to automatically map integration sites to target genes. 

# Getting started
1. Make sure the packages biomaRt, BSgenome* of choice, and [CIMPL](https://github.com/NKI-CCB/cimpl/) v1.0 are installed.
2. Execute the sample code in file 'prog.r' in directory 'r' to check if everything works.
3. Adapt code to use your own input files and configurations.

# Code documentation
#### Directories
- data: sample MuLV and SB transposon data for use with the sample program
- r: KC-RBM scripts, and sample script for running KC-RBM

#### Ouput
KC-RBM returns a data.frame with the following columns:

All fields in the original idata, and additional fields indicating for each insertion (or insertion cluster):
- $transid: to which transcript ID(s) this insertion maps, including all possible transID with all possible mechanisms
- $d2tss: distance to TSS
- $d2tts: distance to TTS
- $dfrac_trans: fractional distance w.r.t. transcript length, relative to TSS (only if insertion is within a transcript).
- $d2css: distance to coding start
- $d2cts: distance to coding end
- $dfrac_cds: fractional distance w.r.t. coding region length, relative to coding start (only if insertion is within a coding region).
- $ensid: to which gene ID(s) the transcript ID(s) belong
- $d2gss: distance to gene start
- $d2gts: distance to gene end
- $dfrac_gene: fractional distance w.r.t. gene length, relative to gene start (only if insertion is within a gene).
- $distance: the distance between the insertion and the TSS of the mapped transcripts(bp, not absolute value) (distance = insertion_locus - TSS)
- $mechanism: by which mechanism the insertion is mapped to a transcript
  - us: upstream-sense (insertion is upstream of a target gene, in sense orientation relative to the gene)
  - ua: upstream-antisense
  - ds: downstream-sense
  - da: downstream-antisense
  - utr3_intron_s: 3'UTR-intron-sense (insertion is in a 3'UTR intron, in sense orientation relative to the gene)
  - utr3_intron_a: 3'UTR-intron-antisense
  - utr5_intron_s: 5'UTR-intron-sense
  - utr5_intron_a: 5'UTR-intron-antisense
  - utr3_exon_s: 3'UTR-exon-sense (insertion is in a 3'UTR exon, in sense orientation relative to the gene)
  - utr3_exon_a: 3'UTR-exon-antisense
  - utr5_exon_s: 5'UTR-exon-sense
  - utr5_exon_a: 5'UTR-exon-antisense
  - other_within_s: any other insertions, within genes and in sense orientation relative to the gene
  - other_within_a: any other insertions, within genes and in antisense orientation relative to the gene
- $gene_mechanism: orientation/location relative to the target gene
  - us: upstream-sense (insertion is upstream of a target gene, in sense orientation relative to the gene)
  - ua: upstream-antisense
  - ds: downstream-sense
  - da: downstream-antisense
  - is: inside-sense (insertion within a gene, in sense orientation relative to the gene)
  - ia: inside-antisense
- $id: numeric identifier of the insertion
- $clusterpeak: location of the nearest KC peak
- $clusterorientation: overall orientation of the cluster consisting of all insertions that have $clusterpeak as the nearest peak
- $ids: $id's of all insertions in this cluster
- $n_insertions: no. of insertions in this cluster

# Reference
De Jong et al. Computational identification of insertional mutagenesis targets for cancer gene discovery. Nucleic Acids Res. 2011 Aug;39(15):e105. [doi: 10.1093/nar/gkr447](http://dx.doi.org/10.1093/nar/gkr447).
