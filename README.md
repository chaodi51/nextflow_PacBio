# Nextflow pipeline for PacBio data
A Nextflow pipeline using PacBio long-read sequencing data

## Purpose
* Vector and plasmid characterization

### Main components

* Basic QC on input reads
* Alignment to a user-supplied reference, a human reference (T2T), and helper and rep-cap plasmids with pbmm2
* Methylation quantification at CpG sites with primrose
* Call variants (lofreq, pbsv) and produce variation table for each reference locus
* Truncation hotspot (read start and end positions) and read length figures
* Chimeric reads detection using Alvis
* Generate positional sequence quality metrics
* Figures, tables, and downloads provided in report on RStudio Connect


## Input
* Unaligned bam of ccs reads
* Reference FASTA
* Various optional control and metadata parameters
