# WheresWalker
WheresWalker is the software pipeline that helps identify regions of homozygosity to assist in experiments that help relate genotype to phenotype in genes with recessive mutations.

# Installation
Whereswalker is written in Bash and Perl, it is cross-platform compatible and no compilation is needed.  To install, run:
```
git clone https://github.com/alekseyzimin/WheresWalker
```
The main executable is whereswalker.sh. Whereswalker also includes as slightly modified for compatibility version of ANNOVAR software.

# Usage
```
whereswalker.sh [arguments]
-m PATH  mutant vcf file MANDATORY
-w PATH  wild type vcf file MANDATORY
-a PATH  annotation gtf file MANDATORY
-g PATH  genome fasta file MANDATORY
-f PATH  GVF file (see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md) optional
-t FLOAT starting threshold, the program will iterate down from this threshold, until an interval is found: default:4.0
-v       verbose switch
-h       this message
```
The GVF (Genome Variation Format) file for the organism is used to optionally screen for known mutations, using Ensembl data.  The GVF format is described here: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gvf.md. From the link:
```
ENSEMBL is providing GVF files for their sequence_alteration data sets at:

ftp://ftp.ensembl.org/pub/current_variation/gvf/
Data from NCBI
The dbVar database at NCBI is providing GVF files for their structural variant data at:

ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/NCBI36/gvf/
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh37/gvf/
ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/gvf/
```
WheresWalker expectc VCF format that is produced by freebayes software with column 10 in format 
