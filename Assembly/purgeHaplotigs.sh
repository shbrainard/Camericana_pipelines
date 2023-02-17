#!/bin/bash
# *** Purge Haplotigs ***
MM2="$UNDERLYINGPATH/minimap2/minimap2.v2.17/minimap2"
SAMTOOLS="$UNDERLYINGPATH/samtools/v1.9/samtools"
PURGEHAPLOTIGS="$UNDERLYINGPATH/.conda/envs/purge_haplotigs/bin/purge_haplotigs"

# STEP 1:
# Align Reads
$MM2 -t 80 -L --MD -ax map-ont $GENOME $READS > $SAM ;

# Convert Sam to Bam, Sort and Index:
$SAMTOOLS sort -m 12G -@ $THREADSby2 -o $SORTED_BAM $SAM ;
$SAMTOOLS index -@ 10 $SORTED_BAM ;

# Generate Histogram
$PURGEHAPLOTIGS hist -bam $SORTED_BAM -genome $GENOME -threads $THREADS ;

echo "[STDOUT]: Exiting 0 now so you can view your histogram output and decide upon cutoffs..."
exit 0;  

# STEP 2:
LOWCUT="3"
MIDCUT="22"
HIGHCUT="190"
# gerneate covdeage stats
$PURGEHAPLOTIGS cov -i ${SORTED_BAM}.gencov -low $LOWCUT -mid $MIDCUT -high $HIGHCUT -o coverage_stats.csv ;
# use coverage stats to purge haplotigs
$PURGEHAPLOTIGS purge -g $GENOME -c coverage_stats.csv -threads 40 ;