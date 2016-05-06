#!/bin/bash

### Description: this is the bash script for determining solo LTRs for zebra finch (taeGut)

### cleaned on 5/5/2016, 
### Yanzhu Ji (yanzhuji20@gmail.com)

## blast against whole genome with 5' and 3' ltr respectively


module load blast
module load BEDTools 

cd $PBS_O_WORKDIR

#### re-run 6/30/15, after removing the "singleton"-ed cluster.

fasta_base_1=taeGut-all_no-CR1-RefSeq_LTRdigest-id80
base=clustered-messy-flanks


### link 5' and 3' ltr fasta files to the current directory.  
ln -s ../harvest-id80/${fasta_base_1}_3ltr.fas .
ln -s ../harvest-id80/${fasta_base_1}_5ltr.fas .
ln -s ../harvest-id80/usearch8-cluster/clustered-aln/clustered-messy-flanks_original.prefa .

### find overlapping between prelim-ERVs and 5'/3' LTRs.
perl ~/perl/fetchFasta.pl ${fasta_base_1}_3ltr.fas ${base}_original.prefa > ${base}_3ltr.fasta
perl ~/perl/fetchFasta.pl ${fasta_base_1}_5ltr.fas ${base}_original.prefa > ${base}_5ltr.fasta

blastn -task blastn \
	-db $RCAC_SCRATCH/taeGut/taeGut-all_noDesc.fna \
	-query ${base}_5ltr.fasta \
	-evalue 0.01 \
	-outfmt 6 \
	-out ${base}_5ltr.fmt6


blastn -task blastn \
	-db $RCAC_SCRATCH/taeGut/taeGut-all_noDesc.fna \
	-query ${base}_3ltr.fasta \
	-evalue 0.01 \
	-outfmt 6 \
	-out ${base}_3ltr.fmt6

perl ~/perl/hitLengthFilter.pl \
	-s ${base}_5ltr.fasta \
	-f ${base}_5ltr.fmt6 \
	-h 1 -c 0.95 \
	-p ${base}_5ltr_hlf.fmt6 \
	> ${base}_5ltr_hlf.txt


perl ~/perl/hitLengthFilter.pl \
	-s ${base}_3ltr.fasta \
	-f ${base}_3ltr.fmt6 \
	-h 1 -c 0.95 \
	-p ${base}_3ltr_hlf.fmt6 \
	> ${base}_3ltr_hlf.txt

perl ~/perl/format-transfer.pl -i fmt6 ${base}_5ltr_hlf.fmt6 
perl ~/perl/format-transfer.pl -i fmt6 ${base}_3ltr_hlf.fmt6


### detect LTRs as long as overlapped by both 3' and 5' ltrs.
cat ${base}_5ltr_hlf.merged.bed ${base}_3ltr_hlf.merged.bed | \
	bedtools sort -i - | bedtools merge -c 4 -o count -i - | awk '$4 > 1' > ${base}_final-ltr.bed
#  2205 hits; the number remained the same for the re-run 6/30/15.


perl ~/perl/column-transfer.pl -i ${base}_final-ltr.bed,1,1 \
	-r $RCAC_SCRATCH/taeGut/chr2acc_chrAdded.txt,2,1 \
	> ${base}_final-ltr_chrID.bed
# 756 lines are skipple; 1449 lines remained.


# remove those overlapped with full-length elements
bedtools intersect -a ${base}_final-ltr_chrID.bed -b ../blastn/${base}_taeGut-all-noDesc_E5_hlf_merged_chrID.bed -f 0.5 -v > ${base}_pure-solo-ltr_chrID.bed
# 1390 lines remained.


