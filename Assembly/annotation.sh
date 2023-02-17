#!/bin/bash

#INPUTS
R1=$1  # these both have to be in .fastq format
R2=$2  # these both have to be in .fastq format
SOFT_MASKED_GENOME="/mnt/projects/Winkler.final.softmasked.fa"
GENOMEF=$(basename $SOFT_MASKED_GENOME)
GENOMEDIR=$(dirname $SOFT_MASKED_GENOME)

# programs and builtins
CONDAENVS="/home/ec2-user/miniconda3/envs"
STAR="$CONDAENVS/read_aligners/bin/STAR"
GSNAP="$CONDAENVS/gmap/bin/gsnap"
GMAP="$CONDAENVS/gmap/bin/gmap"
GMAP_BUILD="$CONDAENVS/gmap/bin/gmap_build"
SAMTOOLS="$CONDAENVS/samtools/bin/samtools" # v1.9
BAMTOOLS="$CONDAENVS/bamtools/bin/bamtools"
OPERATION_N_FILTER="/mnt/projects/genome_assembly/modules/genome_annotation/operation_N_filter.txt"
WORKINGDIR=$(pwd)
DOCKER_RUN_AUGUSTUS="docker run -v $WORKINGDIR:/opt/temp/ --rm sandmansdownfall/augustus:latest"
EXTRINSIC_CONFIG="/root/augustus/config/extrinsic/extrinsic.M.RM.E.W.cfg" # only exists in docker container
SPECIES="tomato"
CPU=$(nproc)
HALFCPU=$(($CPU/2))

*** AUGUSTUS STAGE 1 ***
#cmds 
# build idx
$STAR --runThreadN $HALFCPU --runMode genomeGenerate --genomeDir $WORKINGDIR/genome_masked --genomeFastaFiles $SOFT_MASKED_GENOME ;

# run STAR alignment
$STAR --outSAMattributes NH HI AS nM NM MD jM jI MC --outFilterMismatchNmax 10 --outFilterMatchNminOverLread 0.66 --runThreadN $CPU --genomeDir $WORKINGDIR/genome_masked --readFilesIn $R1 $R2 --readFilesCommand zcat --outSAMtype SAM  ;
mv Aligned.out.sam rnaseq.sam ;

# sort 
$SAMTOOLS view -bS -o rnaseq.bam rnaseq.sam ;
$SAMTOOLS sort -@$CPU -n -o rnaseq.ss.bam rnaseq.bam ;

# augustus annotation run 1
$DOCKER_RUN_AUGUSTUS filterBam --uniq --paired --pairwiseAlignment --in /opt/temp/rnaseq.ss.bam --out /opt/temp/rnaseq.f.bam ;

$SAMTOOLS sort -@$CPU -o rnaseq.fs.bam rnaseq.f.bam ;
$SAMTOOLS view -H rnaseq.fs.bam > header.txt ;

# generate bam hints files showing introns
$DOCKER_RUN_AUGUSTUS bam2hints --intronsonly --in=/opt/temp/rnaseq.fs.bam --out=/opt/temp/hints.gff ;

# generate initial gff
$DOCKER_RUN_AUGUSTUS augustus --species=$SPECIES --extrinsicCfgFile=$EXTRINSIC_CONFIG --alternatives-from-evidence=true --hintsfile=/opt/temp/hints.gff --allow_hinted_splicesites=atac --introns=on --genemodel=complete /opt/temp/$GENOMEF > augustus.gff ;
