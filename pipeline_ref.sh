#!/bin/bash

# New pipeline as of 2013-10-22.

ROOT=/home/quickassist
BWA=$ROOT/genomics/bwa-0.7.8/bwa

REF=$ROOT/genomics/ReferenceMetadata
FASTA=$REF/human_g1k_v37.fasta
FAI=$FASTA.fai
DBSNP=$REF/dbsnp-all.vcf
DBINDEL=$REF/1000G_phase1.indels.b37.vcf
IDIR=$ROOT/genomics/InputFiles
INFILE=HCC1954
ODIR=$ROOT/genomics/OutputFiles
NTHREAD=4
NCTHREAD=4
NBATCH_SIZE=512

HDR=`printf "'@RG\tID:%s\tLB:%s\tSM:%s'" $INFILE $INFILE $INFILE`
#IN1=$IDIR/${INFILE}_1_100reads.fq
IN1=$IDIR/${INFILE}_1_1Mreads.fq
OFILE=$ODIR/$INFILE
SAMFILE=${OFILE}_ref.sam

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
	R=$?
	# Print the end time
	D=`date`
	printf "\n%s\n---- EXIT CODE %s\n---- %s\n%s\n\n\n" $L $R "$D" $L
	if [ $R -ne 0 ]; then
		exit
	fi
}

# Align sequences with BWA
WRAPPER "bwamem" "$BWA mem -t $NTHREAD -Ma -R $HDR $FASTA $IN1 > $SAMFILE"

# DONE!

