#!/bin/bash
DY=$1
WT=$2
THRESH=2
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
set -o pipefail
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
echo "Usage:"
echo "evaluate_mutations.sh [arguments]"
echo "-m <mutant vcf file:string MANDATORY>"
echo "-w <wild type vcf file:string MANDATORY>"
echo "-a <annotation gtf file:string MANDATORY>"
echo "-g <genome fasta file:string MANDATORY>"
echo "-f <GVF file:string optional>"
echo "-t <threshold:float default:2>"
echo "-v verbose switch"
echo "-h help message"
echo ""
}

#parsing arguments
if [[ $# -eq 0 ]];then
usage
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -m|--mutant)
            export DY="$2";
            shift
            ;;
        -w|--wildtype)
            export WT="$2";
            shift
            ;;
        -f|--gvf)
            export GVF_FILE="$2";
            shift
            ;;
        -a|--annotation)
            export ANNOTATION_GTF="$2";
            shift
            ;;
        -g|--genome)
            export GENOME="$2";
            shift
            ;;
        -t|--threshold)
            export THRESH="$2";
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [ ! -s $DY ];then 
  error_exit "Mutant VCF file $DY not found or size zero"
fi
if [ ! -s $WT ];then
  error_exit "Wild Type VCF file $WT not found or size zero"
fi
if [ ! -s $GENOME ];then
  error_exit "Genome sequence FASTA file $GENOME not found or size zero"
fi
if [ ! -s $ANNOTATION_GTF ];then
  error_exit "Annotation GTF file $ANNOTATION_GTF not found or size zero"
fi

DY_VCF=`basename $DY`
WT_VCF=`basename $WT`

if [ ! -e frequencies.success ];then
  log "Computing intervals homozygous in the mutant"
  grep --color=auto -v '^#' $DY | awk '{split($10,a,":");if(a[4]>2 && a[6]>2) print $1" "$2" "$4" "$5" "a[4]" "a[6]}' | $MYPATH/compute_snp_freq.pl > $DY_VCF.DY.txt.tmp &
  PID1=$!;
  grep --color=auto -v '^#' $WT | awk '{split($10,a,":");if(a[4]>2 && a[6]>2) print $1" "$2" "$4" "$5" "a[4]" "a[6]}' | $MYPATH/compute_snp_freq.pl > $WT_VCF.WT.txt.tmp & 
  PID2=$!;
  wait $PID1 $PID2 && mv $DY_VCF.DY.txt.tmp $DY_VCF.DY.txt && mv $WT_VCF.WT.txt.tmp $WT_VCF.WT.txt && \
  touch frequencies.success
fi

if [ ! -e intervals.success ];then
  paste $DY_VCF.DY.txt $WT_VCF.WT.txt |awk '{if(NF==6) print $1" "$2" "$3" "$6}'| $MYPATH/compute_intervals.pl $THRESH > intervals.txt.tmp && \
  mv intervals.txt.tmp intervals.txt && \
  NUM_INTERVALS=`wc -l intervals.txt | awk '{print $1-1}'` && \
  if [ $NUM_INTERVALS -lt 1 ];then
    error_exit "No intervals found, try re-running with a lower threshold (-t), current threshold is $THRESH"
  else
    log "Found $NUM_INTERVALS intervals in intervals.txt"
  fi && \
  touch intervals.success || error_exit "Computing intervals failed"
fi

if [ ! -r prepare.success ];then
  log "Preparing ANNOVAR database"
  $MYPATH/ANNOVAR/tools/gtfToGenePred -genePredExt <(grep -v unknown_transcript $ANNOTATION_GTF) DR_refGene.txt.tmp && \
  mv DR_refGene.txt.tmp DR_refGene.txt && \
  perl $MYPATH/ANNOVAR/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile $GENOME DR_refGene.txt --out DR_refGeneMrna.fa.tmp 2>/dev/null && \
  mv DR_refGeneMrna.fa.tmp DR_refGeneMrna.fa && \
  touch prepare.success || error_exit "Preparing ANNOVAR database failed"
fi

if [ ! -e examine_genes.success ];then
  log "Examining mutations in genes"
  cat <(grep '^#'  $DY) \
  <(cat \
  <(tail -n +2 intervals.txt |perl -ane '{print "$F[0]\t$F[1]\tstart\n$F[0]\t$F[2]\tend\n"}') \
  <(perl -ane '{$h{$F[0]}=1}END{open(FILE,"'$DY'");while($line=<FILE>){chomp($line);@f=split(/\t/,$line,2);print "$line\n" if(defined($h{$f[0]}));}}' intervals.txt) \
  |sort -k1,1 -k2,2n -S 10% | awk 'BEGIN{flag=0;}{if($3=="start"){flag=1}else if($3=="end"){flag=0};if(flag){split($10,a,":");if($1 ~/^#/|| ((a[4]<2 && a[6]>=2) && (length($4)==length($5)) && (length($4)<=2) && $4!="N" && $5!="N")) print $0}}') | uniq >  $DY_VCF.filtered.vcf && \
  perl $MYPATH/ANNOVAR/annovar/table_annovar.pl $DY_VCF.filtered.vcf ./ --vcfinput --outfile annov -buildver DR --protocol refGene --operation g 1>annovar_out.txt 2>&1 && \
  awk -F '\t' '{split($NF,a,":");if($6~/UTR/ || $9=="stoploss" || $9=="stopgain" || $9=="startloss" || $9=="startgain" || $9=="nonsynonymous SNV" || $9=="nonframeshift substitution") print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"a[4]"\t"a[6]}' annov.DR_multianno.txt | \
  sed 's/ /_/g'|  \
  perl -ane 'BEGIN{
    print "Chr\tCoord\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tmut.ref\tmut.alt\twt.ref\twt.alt\tratio\n";
    open(FILE,"'$WT'");
    while($line=<FILE>){
      chomp($line);
      @f=split(/\t/,$line);
      @ff=split(":",$f[-1]);
      $ratio=$ff[3]/($ff[5]+0.00000001);
      $h{"$f[0] $f[1]"}="$ff[3]\t$ff[5]\t$ratio";
    }
  }
  {
    if(defined($h{"$F[0] $F[1]"})){
      @f=split(/\t/,$h{"$F[0] $F[1]"});
      print join("\t",@F),"\t",$h{"$F[0] $F[1]"},"\n" if($f[2]>1 && $f[2]<4);
    }else{
      print join("\t",@F),"\tNA\tNA\tNA\n"
    }
  }' > genes_to_examine.with_WT_freq.txt.tmp && \
  mv  genes_to_examine.with_WT_freq.txt.tmp  genes_to_examine.with_WT_freq.txt && \
  log "Mutations in genes are annotated in genes_to_examine.with_WT_freq.txt" && touch examine_genes.success || error_exit "Examining mutations failed"
fi

if [ ! -e add_sift.success ];then
  if [ -s $GVF_FILE ];then
    log "Adding SIFT information"
    awk  '{ind=-1;for(i=1;i<NF;i++){if($i=="tolerated" || $i=="tolerated_-_low_confidence") ind=i+1;} if(ind>-1){print $1,$4,$(ind-1),$ind}}' $MYPATH/$GVF_FILE |\
    perl -ane '$h{"$F[0] $F[1]"}=$F[3];END{$h{"Chr Coord"}="SIFT";open(FILE,"genes_to_examine.with_WT_freq.txt");while($line=<FILE>){chomp($line);@f=split(/\t/,$line);print "$line\t";if(defined($h{"$f[0] $f[1]"})){print $h{"$f[0] $f[1]"},"\n";}else{print "NA\n"}}}' > genes_to_examine.with_WT_freq.wSIFT.txt.tmp && \
    mv genes_to_examine.with_WT_freq.wSIFT.txt.tmp genes_to_examine.with_WT_freq.wSIFT.txt && \
    log "Mutations in genes with SIFT information added are in genes_to_examine.with_WT_freq.wSIFT.txt" && touch add_sift.success
  fi
fi

if [ -e examine_genes.success ];then
  log "Candidate mutations are in genes_to_examine.with_WT_freq.txt"
fi

if [ -e add_sift.success ];then
  log "Candidate mutations with SIFT scores are in genes_to_examine.with_WT_freq.wSIFT.txt"
fi

if [ -e examine_indels.success ];then
  log "Computing indel bias and PCR targets"
  bash $MYPATH/get_indels.sh $DY $WT && 
  touch examine_indels.success && \
  log "Indel bias index is in mut.indels.index.txt and the PCR targets are in mut.pcr_targets.txt"
fi

