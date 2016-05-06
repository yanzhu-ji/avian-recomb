#!/bin/bash

module load blast
module load BEDTools

cd $PBS_O_WORKDIR

fasta_base_subset=taeGut-all_prelim-ERVs-2 # re-run 4.11.2015
blastn -task blastn \
	-db $RCAC_SCRATCH/taeGut/taeGut-all_noDesc.fna \
	-query $RCAC_SCRATCH/LTRharvest/taeGut/harvest-id80/${fasta_base_subset}.fasta \
	-evalue 0.00001 \
	-outfmt 6 \
	-out ${fasta_base_subset}_taeGut-all-noDesc_E5.fmt6


perl ~/perl/hitLengthFilter.pl  -s $RCAC_SCRATCH/LTRharvest/taeGut/harvest-id80/${fasta_base_subset}.fasta \
	-f ${fasta_base_subset}_taeGut-all-noDesc_E5.fmt6 \
	-h 1 -c 0.9 -p ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf.fmt6 \
	> ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf.txt

perl ~/perl/format-transfer.pl \
	-i fmt6 ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf.fmt6

module load BEDTools

bedtools sort -i ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf.bed | \
	bedtools merge -i - > ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf_merged.bed

perl ~/perl/column-transfer.pl \
	-i ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf_merged.bed,1,1 \
	-r $RCAC_SCRATCH/taeGut/chr2acc_chrAdded.txt,2,1 \
	> ${fasta_base_subset}_taeGut-all-noDesc_E5_hlf_merged_chrID.bed
# 38 lines are skipped, with 100 remaining...

#######################################################################################################
### updated 5/15, blastn using the LTR-RTs (1) with messy flanking regions and (2) within a cluster. ###
### re-run 6/30, after removing the weird "singleton"-ed cluster

ln -s ../harvest-id80/usearch8-cluster/clustered-aln/clustered-messy-flanks.fasta .

fasta_base_subset_2=clustered-messy-flanks
blastn -task blastn \
	-db $RCAC_SCRATCH/taeGut/taeGut-all_noDesc.fna \
	-query ${fasta_base_subset_2}.fasta \
	-evalue 0.00001 \
	-outfmt 6 \
	-out ${fasta_base_subset_2}_taeGut-all-noDesc_E5.fmt6


perl ~/perl/hitLengthFilter.pl  -s ${fasta_base_subset_2}.fasta \
	-f ${fasta_base_subset_2}_taeGut-all-noDesc_E5.fmt6 \
	-h 1 -c 0.9 -p ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf.fmt6 \
	> ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf.txt

perl ~/perl/format-transfer.pl \
	-i fmt6 ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf.fmt6

module load BEDTools

bedtools sort -i ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf.bed | \
	bedtools merge -i - > ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf_merged.bed

perl ~/perl/column-transfer.pl \
	-i ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf_merged.bed,1,1 \
	-r $RCAC_SCRATCH/taeGut/chr2acc_chrAdded.txt,2,1 \
	> ${fasta_base_subset_2}_taeGut-all-noDesc_E5_hlf_merged_chrID.bed
# 11 lines skipped.

#### end of re-run 6/30/15.

