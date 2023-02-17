#!/bin/bash

#*** TASSEL5 ***

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -GBSSeqToTagDBPlugin -i $WORKINGDIR/tmp/Illumina/ -k "'${V2_KEY}'" -e ${ENZYME} -db $DB -c $REQUIRED_TAG_COUNT -kmerLength ${kmerLength} -minKmerL 20 -mnQS ${mnQS} -mxKmerNum 100000000 -endPlugin -runfork1 ;

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -TagExportToFastqPlugin -db $DB -o $WORKINGDIR/tmp/mergedTagCounts/mergedAll.fq -c $REQUIRED_TAG_COUNT -endPlugin -runfork1 ;

${BOWTIE2}/bowtie2 -p 40 --very-sensitive-local -x ${GENOME_INDEX}/${GENOME_BASE} -U $WORKINGDIR/tmp/mergedTagCounts/mergedAll.fq -S $WORKINGDIR/tmp/topm/alignment.sam ;

${TASSEL5}/run_pipeline.pl $MEM100G -fork1 -SAMToGBSdbPlugin -i $WORKINGDIR/tmp/topm/alignment.sam -db $DB -minMAPQ 10 -aProp 0.0 -aLen 20 -endPlugin -runfork1 ;

${SAMTOOLS}/samtools view -bh $WORKINGDIR/tmp/topm/alignment.sam | ${SAMTOOLS}/samtools sort -@ 10 - > $WORKINGDIR/tmp/topm/alignment.bam ;

${SAMTOOLS}/samtools flagstat $WORKINGDIR/tmp/topm/alignment.bam > $WORKINGDIR/tmp/topm/alignment.flagstat ;

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -DiscoverySNPCallerPluginV2 -db $DB -mnMAF 0.01 -mnLCov 0.1 -endPlugin -runfork1 ;

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -SNPQualityProfilerPlugin -db $DB -statFile $WORKINGDIR/tmp/snpCaller/stats.txt -endPlugin -runfork1 ;

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -GetTagSequenceFromDBPlugin -db $DB -o $WORKINGDIR/VCF/retagged/reference_tags.txt -endPlugin -runfork1 ;

${TASSEL5}/run_pipeline.pl $MEMby2 -fork1 -ProductionSNPCallerPluginV2 -i $WORKINGDIR/tmp/Illumina -k "${V2_KEY}" -e ${ENZYME} -kmerLength ${kmerLength} -mnQS ${mnQS} -db $DB -o $WORKINGDIR/VCF/SNPs.mergedAll.vcf -endPlugin -runfork1
# fix a minor error in Tassel5 output files
sed -i 's/##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">/##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n##INFO=<ID=QualityScore,Number=1,Type=Float,Description=\"Quality Score\">/g' $WORKINGDIR/VCF/SNPs.mergedAll.vcf

# generate MAF>0.05 filtered SNP file
$VCFTOOLS --vcf $WORKINGDIR/VCF/SNPs.mergedAll.vcf --remove-indels --maf 0.05 --recode --out $WORKINGDIR/VCF/Filtered/SNPs.mergedAll.MAFgt0.05 ;