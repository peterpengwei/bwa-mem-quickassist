#!/bin/bash

# New pipeline as of 2013-10-22.

# gene1 data path
ROOT=/home/hadoopmaster/genomics

BWA=$HOME/genomics/bwa-mem-quickassist/bwa-0.7.8/bwa
BWA_GOLDEN=$HOME/genomics/bwa_golden

REF=$ROOT/ReferenceMetadata
FASTA=$REF/human_g1k_v37.fasta
FAI=$FASTA.fai
DBSNP=$REF/dbsnp-all.vcf
DBINDEL=$REF/1000G_phase1.indels.b37.vcf
IDIR=$ROOT/InputFiles
INFILE=HCC1954
ODIR=/tmp/cody
NTHREAD=24 # gene1 has 24 threads
NCTHREAD=4
NBATCH_SIZE=300

HDR=`printf "'@RG\tID:%s\tLB:%s\tSM:%s'" $INFILE $INFILE $INFILE`
READ_SIZE=100M
IN1=$IDIR/${INFILE}_1_${READ_SIZE}reads.fq
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
	echo "---------------------" >> flow.log
	echo "A new task start at $D" >> flow.log
	echo "$2" >> flow.log
#	# Run the command using strace
#	eval "strace -T -ttt -o $OFILE.strace.$1 -f -ff $2"
	# Run the command
	eval "$2"
	R=$?
	# Print the end time
	D=`date`
	printf "\n%s\n---- EXIT CODE %s\n---- %s\n%s\n\n\n" $L $R "$D" $L
	echo "End at $D with EXIT CODE $R" >> flow.log
	echo "---------------------" >> flow.log

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
USE_BATCH=1

rm $SAMFILE
if [ $USE_BATCH == 1 ]; then
	WRAPPER "bwamem" "$BWA mem -t $NTHREAD -b $NBATCH_SIZE -Ma -R $HDR $FASTA $IN1 > $SAMFILE"
else
	WRAPPER "bwamem" "$BWA_GOLDEN mem -t $NTHREAD -Ma -R $HDR $FASTA $IN1 > $SAMFILE"
fi

#echo Diff: | tee -a flow.log
#diff $SAMFILE ${ODIR}/_${READ_SIZE}HCC1954_dut.sam | tee -a flow.log

