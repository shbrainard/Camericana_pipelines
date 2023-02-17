#!/bin/bash
#*** RagTag ***

QUERYASSEMBLY="$1"
REFGENOME="$2"  #UNDERLYINGPATH/genomes/Corylus_avellana/GCA_901000735.2_Tombul_CavTom2PMs-1.0/genome.fa
INPUTREADS="$3" # must be in fasta format: 
DATATYPE="ont"
#  -T STR                read type. 'sr', 'ont' and 'corr' accepted for
#                        Illumina, nanopore and error corrected long-reads,
#                        respectively [null]
QUERYBASE=$(basename $QUERYASSEMBLY .fa)
REFBASE=$(basename $REFGENOME .fa)
WORKINGDIR=$(pwd)
UNDERLYINGPATH=""
DIROUT="$WORKINGDIR/${QUERYBASE}_x_${REFBASE}"
SCAFFOLD_DIROUT="$WORKINGDIR/${QUERYBASE}_x_${REFBASE}/ragtag_scaffolded"
RAGTAG="$UNDERLYINGPATH/anaconda3/envs/ragTag_2.0.1/bin/ragtag.py"
NUCMER="$UNDERLYINGPATH/mummer/4.0.0/nucmer"
THREADS=40
SKIPFILE=""
#Skipfile is just a text file with chromosomes to skip
# LR899423.1_pchr01_HIC_scaffold_1
# LR899427.1_pchr05_HIC_scaffold_2
# LR899428.1_pchr06_HIC_scaffold_3

# check for skip file
if [ "$SKIPFILE" != "" ];
then
    echo "[STDOUT]: skipping chromsome correction on:";
    cat $SKIPFILE ;
    CORRECT_OPTIONS=" -j $SKIPFILE "
else
    echo "[STDOUT]: correcting all chromosomes";
    CORRECT_OPTIONS=""
fi
# check for reads 
if [ "$INPUTREADS" == "" ];
then
    echo "[STDOUT]: input reads not found";
    CORRECT_OPTIONS="${CORRECT_OPTIONS}"
else
    echo "[STDOUT]: input reads found";
    CORRECT_OPTIONS="${CORRECT_OPTIONS} -R ${INPUTREADS}"
fi
# run ragtag correct: 
$RAGTAG correct -T ${DATATYPE} ${CORRECT_OPTIONS} -t $THREADS -o $DIROUT $REFGENOME $QUERYASSEMBLY ;


echo "[STDOUT]: entering $DIROUT to perform scaffolding on $DIROUT/ragtag.correct.fasta"
cd $DIROUT ;

# run ragtag scaffold
$RAGTAG scaffold -o $SCAFFOLD_DIROUT -t $THREADS $REFGENOME $DIROUT/ragtag.correct.fasta

# run nucmer to compare the original Refbase assembly against query base:
$NUCMER -p ${QUERYBASE}_x_${REFBASE} -t $THREADS --maxmatch -l 100 -c 500 ${REFGENOME} "$SCAFFOLD_DIROUT/ragtag.scaffold.fasta"

# load delta files into http://assemblytics.com/ 