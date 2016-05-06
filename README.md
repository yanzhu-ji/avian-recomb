This is the readme file for dataset avian_recomb.  Files in this dataset include:
(1) folder R_data_scripts/
    chicken (galGal4) and zebra finch (taeGut) files that are used as input in R, as well as R scripts;
(2) bash scripts (*.sh, with zebra finch as the example) and stringent LTR-RTs in fasta format for both species;
(3) folder perl/
    perl scripts that are called from the bash scripts.

Detailed description:

(1) Files in R_data_scripts/ include:

a. recombination rate files after filtering of dislocated markers
galGal4_new_5-col.out
taeGut_new_5-col_chrID.out

b. BED files of locations of full-length LTR-RTs and solo LTRs
galGal4_relaxed_LTR-RT_chrID.bed
galGal4_relaxed_solo-ltr_chrID.bed

taeGut_relaxed_LTR-RT_chrID.bed
taeGut_relaxed_solo-ltr_chrID.bed

c. BED file with locations of centromeres in chicken
galGal4_centromere.bed

d. R scripts: 
galGal4_recomb_vs_f-s-ratio_upload.R
taeGut_recomb_vs_f-s-ratio_upload.R


(2) bash files and fasta:

To identify LTR-RTs and get the stringent set in the first place:
taeGut-all_LTRharvest-id80_stringent_clean_050516.sh

To get the relaxed set of LTR-RTs:
taeGut_blastn_relaxed_clean_050516.sh

To get solo LTRs:
taeGut_solo-blastn_clean_050516.sh


Stringent set of LTR-RTs in chicken and zebra finch:
galGal4_stringent_LTR-RT.fasta
taeGut_stringent_LTR-RT.fasta


(3) Files in perl/ folder include perl scripts that are used to do a variety of jobs.
Please see description within each file.




