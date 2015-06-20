#!/bin/bash

# New pipeline as of 2013-10-22.

ROOT=/cdsc_nfs/cdsc0/software/spark
BWA=/curr/haoyc/bwa-mem-quickassist/bwa-0.7.8/bwa
#BWA=$ROOT/bwa-0.7.8/bwa
#SAMTOOLS=$ROOT/samtools-0.1.19/samtools
#PICARD=$ROOT/picard-tools-1.79
#GATK=$ROOT/gatk-protected-master/dist/GenomeAnalysisTK.jar

REF=$ROOT/genomics_data/ReferenceMetadata
FASTA=$REF/human_g1k_v37.fasta
FAI=$FASTA.fai
DBSNP=$REF/dbsnp-all.vcf
DBINDEL=$REF/1000G_phase1.indels.b37.vcf
IDIR=$ROOT/genomics_data
INFILE=HCC1954
# ODIR=$ROOT/genomics_data
ODIR=/curr/haoyc
NTHREAD=4
NCTHREAD=4
NBATCH_SIZE=10

HDR=`printf "'@RG\tID:%s\tLB:%s\tSM:%s'" $INFILE $INFILE $INFILE`
#IN1=$IDIR/${INFILE}_1_100reads.fq
IN1=$IDIR/${INFILE}_1_10Mreads.fq
IN2=$IDIR/${INFILE}_2.fq
OFILE=$ODIR/$INFILE
SAMFILE=${OFILE}_dut.sam
UNSORTED=$OFILE.unsorted.bam
BAMFILE=$OFILE.bam
DUPBAM=$OFILE.dup.bam
DUPMETRICS=$OFILE.dup.metrics
TMP=$ROOT/genomics_data
REALN=$OFILE.realn.intervals
REALNBAM=$OFILE.realn.bam
RECAL=$OFILE.recal.grp
FINALBAM=$OFILE.final.bam
FLAGSTAT=$OFILE.flagstat

# Helper function:
#   $1 : tag
#   $2 : command
WRAPPER() {
	L="--------------------"
	L="$L$L$L$L"
	# Print the command and start time
	D=`date`
	printf "%s\n---- %s\n---- %s\n%s\n\n" $L "$D" "$2" $L
#	# Run the command using strace
#	eval "strace -T -ttt -o $OFILE.strace.$1 -f -ff $2"
	# Run the command
	eval "$2"
	echo "$2"
	R=$?
	# Print the end time
	D=`date`
	printf "\n%s\n---- EXIT CODE %s\n---- %s\n%s\n\n\n" $L $R "$D" $L
	if [ $R -ne 0 ]; then
		exit
	fi
}

# One-time preparation of reference genome:
# Create BWA index
#   $BWA index $FASTA
# Create FAI index
#   samtools faidx $FASTA

# Align sequences with BWA
#WRAPPER "bwamem" "$BWA mem -t $NTHREAD -Ma -R $HDR $FASTA $IN1 > $SAMFILE"
#WRAPPER "bwamem" "$BWA --target=ASE mem -t $NTHREAD -b $NBATCH_SIZE -Ma -R $HDR $FASTA $IN1 > $SAMFILE"
WRAPPER "bwamem" "$BWA --target=ASE mem -t $NTHREAD -b $NBATCH_SIZE -Ma -R $HDR $FASTA $IN1"
##WRAPPER "bwamem" "$BWA mem -t $NTHREAD -Ma -R $HDR $FASTA $IN1 $IN2 > $SAMFILE"
##
### Convert SAM to BAM
##WRAPPER "samview" "$SAMTOOLS view -@ $NTHREAD -b -t $FAI -S $SAMFILE -o $UNSORTED"
##
### TEMP: remove SAM file to free disk space
##WRAPPER "" "rm $SAMFILE"
##
### Sort the BAM
##WRAPPER "samsort" "$SAMTOOLS sort -@ $NTHREAD -l 5 $UNSORTED $OFILE"
##
### Index the sorted BAM
##WRAPPER "samindex" "$SAMTOOLS index $BAMFILE"
##
### Mark duplicates in the BAM
##WRAPPER "markdup" "java -Xmx8g -jar $PICARD/MarkDuplicates.jar I=$BAMFILE O=$DUPBAM M=$DUPMETRICS TMP_DIR=$TMP VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true"

### Index the duplicates
##WRAPPER "dupindex" "$SAMTOOLS index $DUPBAM"
##
### Find realignments based on known mutations.
##WRAPPER "realigner" "java -Xmx8g -jar $GATK -T RealignerTargetCreator -nt $NTHREAD -R $FASTA -o $REALN -known:dbsnp,VCF $DBSNP -known:indels,vcf $DBINDEL -I $DUPBAM"
##
### Perform realignments
##WRAPPER "realigner2" "java -Xmx8g -Djava.io.tmpdir=$TMP -jar $GATK -T IndelRealigner -rf NotPrimaryAlignment -R $FASTA -targetIntervals $REALN -known:indels,vcf $DBINDEL -I $DUPBAM -o $REALNBAM"
##
### Index realigned BAM
##WRAPPER "samindex2" "$SAMTOOLS index $REALNBAM"

# Recalibrate quality scores
##WRAPPER "recalib" "java -Xmx4g -jar $GATK -T BaseRecalibrator -I $REALNBAM -R $FASTA -nct $NCTHREAD -rf BadCigar -knownSites:mask,vcf $DBSNP -l INFO --default_platform illumina -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -cov ReadGroupCovariate --disable_indel_quals -o $RECAL"
##
### Print all reads with recalibrated scores into final BAM
##WRAPPER "printreads" "java -Xmx8g -jar $GATK -T PrintReads -R $FASTA -I $REALNBAM -BQSR $RECAL -rf BadCigar -o $FINALBAM"
##
### Index final BAM
##WRAPPER "samindex3" "$SAMTOOLS index $FINALBAM"
##
### Compute final stats
##WRAPPER "flagstat" "$SAMTOOLS flagstat $FINALBAM > $FLAGSTAT"
##
### DONE!
##
