#!/bin/bash
#*** Hi-C ***

# convert Rush-Hazelnut-HiC_S1_L002_R1_001.fastq.gz --> Rush-Hazelnut-HiC
SAMPNAME=$(echo "$FASTQR1_BASE" | sed -E 's/_S.*.fastq.gz//g')
RESTRICTION_ENZ='Arima' 
# RESTRICTION_ENZ could be any of the following
#     'HindIII'     : 'AAGCTT',
#     'DpnII'       : 'GATC',
#     'MboI'        : 'GATC',
#     'Sau3AI'      : 'GATC',
#     'Arima'       : [ 'GATC', 'GANTC' ],

#DIR STRUCTURE
WORKINGDIR=`pwd`
UNDERLYINGPATH=""
GENOME="/path/to/soft/masked/genome.fa"
GENOMEBASE=$(basename $GENOME)
GENOMEDIR=$(dirname $GENOME)
GENOMEID="${GENOMEBASE%.fa}"
CHROMSIZES="$GENOMEDIR/${GENOMEID}.chrom.sizes"
REFFASTA="$WORKINGDIR/references/${GENOMEBASE}"
RESTRICTION_SITE_DIR="$WORKINGDIR/restriction_enz_sites"
CHECK_SITE_DIR=$(ls $RESTRICTION_SITE_DIR/*.txt | wc -l)
JUICER="$UNDERLYINGPATH/progs/juicer/CPU/juicer.sh"
export PATH="$UNDERLYINGPATH/progs/parallel-20181122/src/:$PATH"

#TOOLS:
BWA="$UNDERLYINGPATH/bwa/v0.7.17/bwa"  #0.7.17-r1188
SAMTOOLS="$UNDERLYINGPATH/bin/samtools"
FASTP="$UNDERLYINGPATH/fastp/v0.20.1/fastp"
PY27="$UNDERLYINGPATH/anaconda3/envs/PacbioSuite_py27/bin/python2.7"
JUICER_GENSITEPOS="$UNDERLYINGPATH/progs/juicer/misc/generate_site_positions.py"
THREED_DNA="$UNDERLYINGPATH/progs/3d-dna/run-asm-pipeline.sh"


#setup dir structure 
mkdir -p "$WORKINGDIR/references" "$WORKINGDIR/$SAMPNAME/fastq" "$RESTRICTION_SITE_DIR" ;
# use fastp to QC and place your new data into correct location
$FASTP --thread 16 --html $WORKINGDIR/$SAMPNAME/fastq/${SAMPNAME}.fastp.QC.html --json $WORKINGDIR/$SAMPNAME/fastq/${SAMPNAME}.fastp.QC.json --in1 $FASTQR1 --in2 $FASTQR2 --out1 $WORKINGDIR/$SAMPNAME/fastq/${SAMPNAME}_R1.fastq.gz --out2 $WORKINGDIR/$SAMPNAME/fastq/${SAMPNAME}_R2.fastq.gz ; 
rm $WORKINGDIR/$SAMPNAME/fastq/${SAMPNAME}.fastp.QC.json ;

# setup scripts
ln -s $GENOME $REFFASTA
cd $WORKINGDIR/references/ ;
$SAMTOOLS faidx $GENOMEBASE ;
$BWA index $GENOMEBASE ;
cd $WORKINGDIR ;


$PY27 $JUICER_GENSITEPOS $RESTRICTION_ENZ ${GENOMEID} $REFFASTA ;
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${GENOMEID}_${RESTRICTION_ENZ}.txt > $CHROMSIZES
cd $WORKINGIDIR ;


# Running Juicer:
rm -r $WORKINGDIR/$SAMPNAME/aligned/ $WORKINGDIR/$SAMPNAME/splits/ $WORKINGDIR/$SAMPNAME/HIC_tmp/ ;
echo "[STDOUT]: $JUICER -p $CHROMSIZES -y $WORKINGDIR/restriction_enz_sites/${GENOMEID}_${RESTRICTION_ENZ}.txt -z $REFFASTA -d $WORKINGDIR/$SAMPNAME -D $WORKINGDIR -t 80 " ;
$JUICER -p $CHROMSIZES -y $WORKINGDIR/restriction_enz_sites/${GENOMEID}_${RESTRICTION_ENZ}.txt -z $REFFASTA -d $WORKINGDIR/$SAMPNAME -D $WORKINGDIR -t 80 ;

# Run 3D DNA
$THREED_DNA --mode 'haploid' --input 15000 --rounds 3 $REFFASTA $WORKINGDIR/$SAMPNAME/aligned/merged_nodups.txt ;
