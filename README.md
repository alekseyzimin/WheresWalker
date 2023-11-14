# WheresWally
WheresWally is the software pipeline that helps identify regions of homozygosity to assist in experiments that help relate genotype to phenotype in genes with recessive mutations.

# Installation

# Usage
```
evaluate_mutations.sh [arguments]
-m <mutant vcf file>:path MANDATORY
-w <wild type vcf file>:path MANDATORY
-a <annotation gtf file>:path MANDATORY
-g <genome fasta file>:path MANDATORY
-f <GVF file:string optional>
-t <starting threshold, the program will iterate down from this threshold, until an interval is found>:float default:4
-v verbose switch
-h help message
```
