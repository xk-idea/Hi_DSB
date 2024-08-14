#!/bin/bash 

# Run this via the following: 
# qsub -v BEDFILE=bedfile -v BIGWIGS=bigwig.list -v OUT=out.txt score_bedfile_with_bigwigs.sh 

#$ -N bigwig 
#$ -cwd
#$ -j y 
#$ -l h_rt=1:00:00
#$ -l h_vmem=4G
# $ -pe sharedmem 1 

scratch=~/scratch
mkdir -p $scratch/$JOB_ID
tmpdir=$scratch/$JOB_ID

awk '{x++; print $1"\t"$2"\t"$3"\t"x}' $BEDFILE | sort -k1,1 -k2,2n > $tmpdir/bedfile.bed 
cp $tmpdir/bedfile.bed $tmpdir/fout.txt
echo -e "chr\tstart\tend\tid" > $tmpdir/header.txt 

for bigwig in `cat $BIGWIGS`;
do 
	label=`basename $bigwig .bigWig`
	label=`basename $bigwig .bw`
	echo "bigWigAverageOverBed $bigwig $tmpdir/bedfile.bed $tmpdir/out.tab"
	bigWigAverageOverBed $bigwig $tmpdir/bedfile.bed $tmpdir/out.tab
	cut -f5 $tmpdir/out.tab | paste $tmpdir/fout.txt - > $tmpdir/tmp 
	mv $tmpdir/tmp $tmpdir/fout.txt
	echo $label | paste $tmpdir/header.txt - > $tmpdir/htmp 
	mv $tmpdir/htmp $tmpdir/header.txt  
done 
cat $tmpdir/header.txt $tmpdir/fout.txt > $OUT

