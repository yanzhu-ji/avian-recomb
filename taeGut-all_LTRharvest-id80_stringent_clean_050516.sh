#!/bin/bash

# run LTRharvest and LTRdigest on galGal4 assembled chromosomes

module load HMMER
module load EMBOSS

cd $PBS_O_WORKDIR

species=taeGut
### 1.1 LTRharvest
gt suffixerator -indexname ${species} -db $RCAC_SCRATCH/${species}/${species}-all_noDesc.fna -tis -suf -lcp -des -ssp -sds -dna 


gt ltrharvest \
	-out ${species}-all_LTRharvest-id80.fasta \
	-gff3 ${species}-all_LTRharvest-id80.gff3 \
	-similar 80 \
	-md5 yes \
	-seqids yes \
	-index ${species}

gt gff3 -sort ${species}-all_LTRharvest-id80.gff3 > ${species}-all_LTRharvest-id80_sorted.gff3
gt md5_to_id ${species}-all_LTRharvest-id80_sorted.gff3 > ${species}-all_LTRharvest-id80_sorted_ids.gff3

### 1.2 CR1 filter

### 2.1 combine LTRharvest and CR1 filter results.  See CR1-RefSeq_filter/
### galGal4-all_LTRharvest-id80_sorted_no-CR1.gff3 is obtained from CR1-RefSeq_filter/ 3.25.2015



### 2.2 LTRdigest on combined LTRharvest and CR1-RefSeq filter results
gt gff3 -sort ${species}-all_LTRharvest-id80_sorted_no-CR1-RefSeq.gff3 > ${species}-all_LTRharvest-id80_no-CR1-RefSeq_sorted.gff3

gt ltrdigest -outfileprefix ${species}-all_no-CR1-RefSeq_LTRdigest-id80 \
	-seqnamelen 50 \
	-uboxdist 15 \
	-uboxlen 3 10 \
	-pptuprob 0.7 \
	-trnas /scratch/lustreC/j/ji20/${species}/${species}-tRNA.fa \
	-seqfile /scratch/lustreC/j/ji20/${species}/${species}-all_noDesc.fna \
	${species}-all_LTRharvest-id80_no-CR1-RefSeq_sorted.gff3 \
	${species}-all > ${species}-all_no-CR1-RefSeq_LTRdigest-id80.gff3

perl ~/perl/filterNN.pl -s ${species}-all_no-CR1-RefSeq_LTRdigest-id80_complete.fas -n 5

### clustering using USEARCH
module load usearch/8.0.1623

base=taeGut-all_no-CR1-RefSeq_LTRdigest-id80_complete_filtered-NN
ln -s ../${base}.fasta .

### cluster all 1452 seqs and save alignment information
usearch8 -cluster_fast ${base}.fasta \
	-id 0.8 \
	-uc ${base}_id80.uc \
	-msaout ${base}_id80-msa.aln

mkdir sgt
for i in ${base}_id80-msa.aln*; do
	number_seqs=`grep -c ">" ${i}`
	if [  "${number_seqs}" -eq "1"  ]; then
	cat $i >> ${base}_sgt.aln  # seemed like putting another tab here would cause problems..not sure why...
	mv $i sgt/
	fi
done

### extend the alignment by 50 bp per side
grep "^C" *.uc |awk '$3 > 2' | cut -d":" -f2 > ${base}_id80_clustered-n3.uc
grep "^C" *.uc |awk '$3 > 2' | cut -f9 > ${base}_id80_clustered-n3-centroids.list

mkdir clustered-aln

# as clip_bed.pl has been modified to print out the name list (-p option), incorporate this step earlier...
for line in `cat ${base}_id80_clustered-n3-centroids.list`; do
	grep $line ${base}_id80.uc | cut -f9 |sort | uniq | sed 's/_\([0-9]\+\)/\t\1/g' > clustered-aln/${line}_cluster-seqs.bed

	perl ~/perl/clip_bed.pl -i clustered-aln/${line}_cluster-seqs.bed \
	-s $RCAC_SCRATCH/taeGut/taeGut-all_noDesc.fna \
	-b 50 \
	-a 50 \
	-p clustered-aln/${line}.txt \
	> clustered-aln/${line}_cluster-seqs_f50.fasta

	usearch8 -cluster_fast clustered-aln/${line}_cluster-seqs_f50.fasta \
	-id 0.75 -uc clustered-aln/${line}_cluster-seqs_f50_uc-id75.uc \
	-msaout clustered-aln/${line}_cluster-seqs_f50_uc-id75-msa.aln

done

cat *.txt > clustered-seqs.txt

##### pfam domains
### link centroids to hmm hits


for line in `cat ${base}_id80_clustered-n3-centroids.list`; do
	grep ${line}_ ../taeGut-all_no-CR1-RefSeq_LTRdigest-id80_NN-not-both_prt.hmmsc.domtblout >> ${base}_id80_clustered-n3-centroids.hmmsc.domtblout
	grep ${line}_ ../taeGut-all_no-CR1-RefSeq_LTRdigest-id80_NN-both_prt.hmmsc.domtblout >> ${base}_id80_clustered-n3-centroids.hmmsc.domtblout
done

perl ~/perl/ERV-filter-v3.pl -info ~/database/info.txt \
	-o ${base}_id80_clustered-n3-centroids_ERV-filter-n3 \
	-p -n 3 -e 0.01 \
	${base}_id80_clustered-n3-centroids.hmmsc.domtblout

cat ${base}_id80_clustered-n3-centroids_ERV-filter-n3.out? > ${base}_id80_clustered-n3-centroids_ERV-filter-n3.out

perl ~/perl/column-transfer.pl -i ${base}_id80_clustered-n3-centroids.list,1,2 \
	-r ${base}_id80_clustered-n3-centroids_ERV-filter-n3.out,2,4 \
	> ${base}_id80_clustered-n3-centroids_pfam.txt

for line in `cat ${base}_id80_clustered-n3-centroids.list`; do
	grep ">" $i clustered-aln/${line}_cluster-seqs_f50.fasta > ${line}_cluster-seqs_f50.prefa
done

### filter out sequences by eyeballing...
### after eyeballing, mark in each prefa file of each cluster the sequences do not qualify (incomplete).  
### Also, write to the file "clustered-messy-flanks_centroids.prefa" the centroids of clusters that have messy flanks.
### In zebra finch, only 5 clusters have messy flanks (totaling 28 seqs)

cd clustered-aln

for line in `cat clustered-messy-flanks_centroids.prefa`; do 
	cat ${line}_clustered-messy-flanks.prefa
done >> clustered-messy-flanks.prefa

# check if the number of sequences match
grep "^>" clustered-messy-flanks.prefa | wc -l


grep "^>" clustered-messy-flanks.prefa > clustered-messy-flanks_seqs.prefa

cat *.txt > clustered-seqs.txt

sed -i 's/^>//' clustered-messy-flanks_seqs.prefa -i 

perl ~/perl/column-transfer.pl -i clustered-messy-flanks_seqs.prefa,1,1 \
	-r clustered-seqs.txt,2,1 \
	> clustered-messy-flanks_original.prefa

perl ~/perl/fetchFasta.pl ../../taeGut-all_no-CR1-RefSeq_LTRdigest-id80_complete_filtered-NN.fasta \
	clustered-messy-flanks_original.prefa > clustered-messy-flanks.fasta