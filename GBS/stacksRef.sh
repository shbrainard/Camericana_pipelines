#!/bin/bash

#*** STACKS 2.2 ***
# simple internal python code that splits combined key files to single key files and outputs them to a directory called $WORKINGDIR/tmp/barcodes
python3 $STACKS_KEYPLITTER --input $KEY_FILE ;
ENZ1=$(echo $ENZYME | sed -E 's/(.+)-(.+)/\1/g' | sed -r 's/\w+/\l&/g')
ENZ2=$(echo $ENZYME | sed -E 's/(.+)-(.+)/\2/g' | sed -r 's/\w+/\l&/g')
#ApeKI --> apeKI
#PstI  --> pstI
if [ "$ENZ2" != "$ENZ1" ]; 
then
	RADENZ="--renz_1 $ENZ1 --renz_2 $ENZ2"
else
	RADENZ="-e $ENZ1"
fi

for BC_file in *BarCode.txt ;
do
    ${STACKS}/process_radtags --barcode-dist-1 0 -t 64 -f $WORKINGDIR/tmp/Illumina/\${BC_file%BarCode.txt}fastq.gz --paired --inline_null --interleaved -b $BC_file -o $WORKINGDIR/tmp/demultiplexed_fastqs $RADENZ -c -q ;
done

# align files:
${BWA}/bwa mem -t 40 $BWA_GENOME_PATH $fastq ${fastq%.1.fq.gz}.2.fq.gz  | ${SAMTOOLS}/samtools view -@ 2 -hb - | ${SAMTOOLS}/samtools sort -@ 4 - > $WORKINGDIR/tmp/aligned/\${fastq%.1.fq.gz}.bam ;

# find alignment rates:
echo "filename total_reads mapping_rate %_mapped" > algn.rates.txt ;
for bam in *.bam;
do 
    TOTAL=$(${SAMTOOLS}/samtools flagstat -@ 40 $bam | head -1 | sed -E 's/^([0-9]+)\s\+.*/\1/g')
    MAPPED=$(${SAMTOOLS}/samtools flagstat -@ 40 $bam | grep 'mapped (' | sed -E 's/([0-9]+).*\(([0-9]+.[0-9]+).*/\1\t\2/g')
    echo $bam $TOTAL $MAPPED >> algn.rates.txt 
done 

# generate popmap with any file > 5Mb reads...
ls -l $WORKINGDIR/tmp/aligned/ | awk -v OFS='\t' '$5>=5000000{print $9,1}' | sed s/.bam//g > $WORKINGDIR/popmap ;

# call gstacks
${STACKS}/gstacks -I $WORKINGDIR/tmp/aligned/ -M $WORKINGDIR/popmap -O $WORKINGDIR/stacks_out/ -t 80

# call populations module
${STACKS}/populations -P $WORKINGDIR/stacks_out/ -M $WORKINGDIR/popmap -r 0.65 --structure --vcf --genepop --fstats --smooth --hwe -t 80

# cleanup
mkdir -p $WORKINGDIR/VCF ;
cd $WORKINGDIR/VCF ;
ln $WORKINGDIR/tmp/stacks_out/populations.snps.vcf ./SNPs.mergedAll.vcf ;
mv $WORKINGDIR/stacks_out $WORKINGDIR/tmp ;


# generate MAF>0.05 filtered SNP file
$VCFTOOLS --vcf $WORKINGDIR/VCF/SNPs.mergedAll.vcf --remove-indels --maf 0.05 --recode --out $WORKINGDIR/VCF/Filtered/SNPs.mergedAll.MAFgt0.05 ;
