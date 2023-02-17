#!/bin/bash

#*** Racon Error Correction ***
#Tool Variables
RACON="$UNDERLYINGPATH/progs/racon/build/bin/racon"
Minimap2="$UNDERLYINGPATH/minimap2/2.17/minimap2"
bwa="$UNDERLYINGPATH/bwa/v0.7.17/bwa"
PBINDEX="$UNDERLYINGPATH/anaconda3/envs/PacbioSuite_py27/bin/pbindex"
BAM2FQ="$UNDERLYINGPATH/anaconda3/envs/PacbioSuite_py27/bin/bam2fastq"

for iteration in {1..3};
do
    echo "[STDOUT]: running iteration # $iteration of Racon Polishing" ;
    if [ $iteration -eq 1 ];
    then
        if [ $isShortRead -eq 0 ];
        then
            # perform alignment, required by Racon
            $Minimap2 --secondary=no -t 40 -x sr $ASSEMBLY $fastq > ${fastq%.fq.gz}.${iteration}.paf ;
            # skipping -L (long cigar format >65535 ops) --MD (include md tag) -c (cigar output in paf), -N (# of secondary alignments to keep default [5]), --secondary=no no secondary alignments accepted, doesn't reduce alignment accuracy
        else
            # perform alignment, required by Racon
            # assumes pacbio CCS is coming in to polish 
            $Minimap2 -N 5 -t 40 -x map-pb $ASSEMBLY $fastq > ${fastq%.fq.gz}.${iteration}.paf ;
            # skipping -L (long cigar format >65535 ops) --MD (include md tag) -c (cigar output in paf), -N (# of secondary alignments to keep default [5]) 
        fi

        # run Racon to polish current assembly
        ASSEMBLYBASE=$(basename $ASSEMBLY)
        $RACON --threads 40 $fastq ${fastq%.fq.gz}.${iteration}.paf $ASSEMBLY > ${ASSEMBLYBASE%.fa}.racon.${iteration}.fa ;
        PREVASSEMBLY="${ASSEMBLYBASE%.fa}.racon.${iteration}.fa" ;

    else
        if [ $isShortRead -eq 0 ];
        then
            # perform alignment, required by Racon
            $Minimap2 --secondary=no -t 40 -x sr $PREVASSEMBLY $fastq > ${fastq%.fq.gz}.${iteration}.paf ;
            # skipping -L (long cigar format >65535 ops) --MD (include md tag) -c (cigar output in paf), -N (# of secondary alignments to keep default [5]), --secondary=no no secondary alignments accepted, doesn't reduce alignment accuracy
        else
            # perform alignment, required by Racon
            $Minimap2 -N 5 -t 40 -x map-pb $PREVASSEMBLY $fastq > ${fastq%.fq.gz}.${iteration}.paf ;
            # skipping -L (long cigar format >65535 ops) --MD (include md tag) -c (cigar output in paf)
        fi
        # run Racon to polish the previous racon assembly
        $RACON --threads 40 $fastq ${fastq%.fq.gz}.${iteration}.paf $PREVASSEMBLY > ${PREVASSEMBLY%.racon.?.fa}.racon.${iteration}.fa ;
        PREVASSEMBLY="${PREVASSEMBLY%.racon.?.fa}.racon.${iteration}.fa" ;
    fi
done