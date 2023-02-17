#!/bin/bash

#*** Assembly ***

READPATH="/path/to/pbconverted/HIFI.reads.fastq"
ESTGENOMESIZE=400
UNDERLYINGPATH=""
PORECHOP="$UNDERLYINGPATH/porechop/0.2.4/porechop"
FILTLONG="$UNDERLYINGPATH/filtlong/0.2.0/filtlong"
minLength=5000


# Pacbio filtering
    $FILTLONG --min_length $minLength $InputFile | gzip > filtered.fq.gz ;

# ONT filtering:
    $PORECHOP --discard_middle --threads 40 --format auto --input $InputFile --output porechopped.fq.gz ;
    # filter by minlength, quality and provide maximum bp coverage...
    $FILTLONG --target_bases $CoverageGoal --min_length $minLength --mean_q_weight 10 porechopped.fq.gz | gzip > filtered.fq.gz ;

HIFIASM="$UNDERLYINGPATH/progs/hifiasm/hifiasm"
$HIFIASM --primary --hg-size $ESTGENOMESIZE -t 160 -o assembly.hifiasm filtered.fq.gz ;

awk '/^S/{print ">"$2;print $3}' assembly.hifiasm.p_ctg.gfa > assembly.hifiasm.fa ;

