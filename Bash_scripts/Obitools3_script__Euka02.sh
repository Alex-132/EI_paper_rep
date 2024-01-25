#!/bin/bash

########################################################################################################
# ObiTools3 Script for Metabarcoding
# original script by Stefanie Knell
# modified by Juliane Romahn (6.4.23: tested without assigning on ObiTool3 v3.0.1b13)
####
# for more details check
# ObiTools3 : https://git.metabarcoding.org/obitools/obitools3
# ObiTools3 Tutorial : https://git.metabarcoding.org/obitools/obitools3/-/wikis/Wolf-tutorial-with-the-OBITools3
##########################################################################################################

########################################## IMPORTANT ##########################################
# an this server pipeline has to be executed in an conda environment, for this extra permissions can be neccessary
# for this execute first: conda activate env_jromahn
##############
# also filenames shoult never contain whitespaces!!!
###########################################################################################

########################################## Dependencies ##########################################
# taxdump.tar.gz has to be in working dir
# obitools3 and if blast assigned wanted blast has to be installed
#########
# EMBL is with the current version of obitools not working (v3.0.1b13)
################################################################################################################

########################################## Parameters ##########################################
# Input
dms="../data/Euka02"    # project name ( no whitespace, underscores instead!!)
read1="../data/AXZS-20230511__230524_A00902_B_L1-2_AXZS-9_R1.fastq.gz" ## path to read1 file
read2="../data/AXZS-20230511__230524_A00902_B_L1-2_AXZS-9_R2.fastq.gz" ## path to read2 file
ngs_file="../data/ngsfilter_GBM_JR_#03_280423__euka02.tsv" # path to ngs file ( has to be tab separated)
#####

##############################
### Cleaning steps
remove_low_alignment="J" #Remove sequences with low overlap alignment score? [J/N] -> has to be capital letters, 
minimum_alignment_score="0.8"  ## alingment score threshold if removing wanted - > recommended 0.8
####
denoise="J" #Denoise the sequence dataset regarding length and count of sequences? [J/N] 
minimum_length="80" # minimum amplicon length
maximum_length="250" # maximum amplicon length
minimum_count="10" # Keep only the sequences having a count greater or equal to X


################################################################################################################

###########################################################################################################
############################## don't change anything after here ############################################
###########################################################################################################

#name project
#read -p "Name of the projekt: " dms
#echo ""
now=$(date)
work_dir=$(pwd)

mkdir ${dms}_results
touch ${dms}_results/read_me.txt

echo "Project name: $dms" > ${dms}_results/read_me.txt
echo "Working directory: $work_dir" >> ${dms}_results/read_me.txt
echo "Forward reads: $read1" >> ${dms}_results/read_me.txt
echo "Reverse reads: $read2" >> ${dms}_results/read_me.txt
echo "NGSfilter: $ngs_file" >> ${dms}_results/read_me.txt
echo "Date: $now" >> ${dms}_results/read_me.txt
echo "ObiTools Version: " >> ${dms}_results/read_me.txt
obi --version >>  ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

#Import the sequencing data into a DMS
obi import $read1 $dms/reads1
obi import $read2 $dms/reads2
obi import --ngsfilter-input $ngs_file $dms/ngsfile


#Recover the full sequences from the partial forward and reverse reads
obi alignpairedend -R $dms/reads2 $dms/reads1 $dms/aligned_reads


#Alignment statisticsSILVA_138.1_LSUParc_tax_silva.fasta
obi stats -a score_norm $dms/aligned_reads > ${dms}_results/01_alignment_statistics.txt
obi stats -m score_norm $dms/aligned_reads >> ${dms}_results/01_alignment_statistics.txt
obi stats -M score_norm $dms/aligned_reads >> ${dms}_results/01_alignment_statistics.txt
obi stats -v score_norm $dms/aligned_reads >> ${dms}_results/01_alignment_statistics.txt
obi stats -s score_norm $dms/aligned_reads >> ${dms}_results/01_alignment_statistics.txt

#clean txt with alignment statistics
sed -i -r 's/score_norm/alignment_score/g' ${dms}_results/01_alignment_statistics.txt 
sed -i -r 's/\scount+\stotal//g' ${dms}_results/01_alignment_statistics.txt 
sed -i -r 's/^(\s)+\s//g' ${dms}_results/01_alignment_statistics.txt
sed -i -r 's/\s\w+\s\w+//g' ${dms}_results/01_alignment_statistics.txt

printf  "01_alignment_statistics.txt \t Overview over alignment scores for the sequences after pairing of forward and reverse reads." >> ${dms}_results/read_me.txt
#echo "Overview over alignment scores for the sequences after pairing of forward and reverse reads." >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

#Remove sequences with low alignment score?
echo ""
cat ${dms}_results/01_alignment_statistics.txt
echo ""
#read -p "Remove sequences with low overlap alignment score? [J/N] " remove_low_alignment

if [ "$remove_low_alignment" == "J" ]
then
    #read -p "Minimum alignment score: " minimum_alignment_score      #(0.8)
    echo ""
    #Remove sequence records with a low overlap alignment score
    obi grep -p "sequence['score_norm'] > $minimum_alignment_score" $dms/aligned_reads $dms/good_sequences
    #Assign each sequence record to the corresponding sample/marker combination
    obi ngsfilter -t $dms/ngsfile -u $dms/unidentified_sequences $dms/good_sequences $dms/identified_sequences
elif [ "$remove_low_alignment" == "N" ]
then
    echo ""
    #Assign each sequence record to the corresponding sample/marker combination
    obi ngsfilter -t $dms/ngsfile -u $dms/unidentified_sequences $dms/aligned_reads $dms/identified_sequences
fi


#Export information about identified and unidentified reads
obi stats -c sample $dms/identified_sequences > ${dms}_results/02_identified_sequences.txt

obi stats -c forward_tag $dms/unidentified_sequences > ${dms}_results/03_unidentified_sequences.txt
echo "" >> ${dms}_results/03_unidentified_sequences.txt
obi stats -c reverse_tag $dms/unidentified_sequences >> ${dms}_results/03_unidentified_sequences.txt

printf "02_identified_sequences.txt \t Number of sequences that could be assigned to their sample" >> ${dms}_results/read_me.txt
#echo "Number of sequences that could be assigned to their sample" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

printf "03_unidentified_sequences.txt \t Sequences that were not assigned to a sample" >> ${dms}_results/read_me.txt
#echo "Sequences that were not assigned to a sample" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt


#Dereplicate reads into unique sequences
obi uniq -m sample $dms/identified_sequences $dms/dereplicated_sequences

#Statistics on length and count of sequences
obi annotate -k COUNT -k MERGED_sample $dms/dereplicated_sequences $dms/cleaned_metadata_sequences

obi stats -a seq_length $dms/dereplicated_sequences  > ${dms}_results/04_length_sequences.txt
sed -i -r 's/^(\s)+\s//g' ${dms}_results/04_length_sequences.txt
obi stats -c seq_length $dms/dereplicated_sequences | sort -n >> ${dms}_results/04_length_sequences.txt
echo ""
cat ${dms}_results/04_length_sequences.txt
echo ""
obi stats -a COUNT $dms/cleaned_metadata_sequences  > ${dms}_results/05_count_sequences.txt
obi stats -c COUNT $dms/cleaned_metadata_sequences | sort -n >> ${dms}_results/05_count_sequences.txt
echo ""
cat ${dms}_results/05_count_sequences.txt | head -20 

printf "04_length_sequences.txt \t Overview over the length of unique sequences after dereplication" >> ${dms}_results/read_me.txt
#echo "Overview over the length of unique sequences after dereplication" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

printf "05_count_sequences.txt \t Overview over the number of occurences of unique sequences after dereplication" >> ${dms}_results/read_me.txt
#echo "Overview over the number of occurences of unique sequences after dereplication" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt


#Denoise the sequence dataset?
echo ""
#read -p "Denoise the sequence dataset regarding length and count of sequences? [J/N] " denoise

if [ "$denoise" == "J" ]
then
    #read -p "Minimum sequence length: " minimum_length                #(80)
    #read -p "Minimum number of Sequences: " minimum_count             #(10)
    echo ""
    obi grep -p "len(sequence)>=$minimum_length and len(sequence)<= $maximum_length and sequence['COUNT']>=$minimum_count" $dms/cleaned_metadata_sequences $dms/denoised_sequences
    obi clean -s MERGED_sample -r 0.05 -H $dms/denoised_sequences $dms/cleaned_sequences
elif [ "$denoise" == "N" ]
then
    obi clean -s MERGED_sample -r 0.05 -H $dms/cleaned_metadata_sequences $dms/cleaned_sequences
fi

##exporting sequences
obi export --fasta-output $dms/cleaned_sequences > ${dms}_results/11_cleaned_sequences.fasta
printf "11_cleaned_sequences.fasta \t Resulting sequences from ObiTools3, for simplicity called cleaned but it isn't if not wanted" >> ${dms}_results/read_me.txt

obi export --tab-output $dms/cleaned_sequences > ${dms}_results/12_cleaned_sequences_table.tsv
printf "12_cleaned_sequences_table.tsv \t Resulting table file from ObiTools3, for simplicity called cleaned but it isn't if not wanted" >> ${dms}_results/read_me.txt

obi export --tab-output $dms/cleaned_sequences > ${dms}_results/13_count_table.tsv
printf "13_count_table.tsv \t Resulting count table from ObiTools3 usable for mothur , for simplicity called cleaned but it isn't if not wanted" >> ${dms}_results/read_me.txt




#Overview Pipeline
obi history -d $dms > ${dms}_results/15_overview_pipeline.dot
dot -Tpng ${dms}_results/15_overview_pipeline.dot > ${dms}_results/16_overview_pipeline.png

obi ls $dms | sed 's/\sDate.*count://g' > ${dms}_results/17_overview_sequences.txt
sed -i '$d' ${dms}_results/17_overview_sequences.txt
sed -i '$d' ${dms}_results/17_overview_sequences.txt

printf "15_overview_pipeline.dot \t Overview graphic for all analysis steps performed with ObiTools in dot format" >> ${dms}_results/read_me.txt
#echo "Overview graphic for all analysis steps performed with ObiTools in dot format" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

printf "16_overview_pipeline.png \t Overview graphic for all analysis steps performed with ObiTools in png format" >> ${dms}_results/read_me.txt
#echo "Overview graphic for all analysis steps performed with ObiTools in png format" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

printf "17_overview_sequences.txt \t Overview over the number of sequences at each step of analysis performed with ObiTools" >> ${dms}_results/read_me.txt
#echo "Overview over the number of sequences at each step of analysis performed with ObiTools" >> ${dms}_results/read_me.txt
echo "" >> ${dms}_results/read_me.txt

echo "" >> ${dms}_results/read_me.txt
echo "More detailed information following:" >> ${dms}_results/read_me.txt
if [ "$remove_low_alignment" == "J" ]
then
    echo "Sequences with an alignment score lower then $minimum_alignment_score were removed." >> ${dms}_results/read_me.txt
fi

if [ "$denoise" == "J" ] 
then
    echo "After dereplicating reads into unique sequences, all sequences shorter than $minimum_length bp , longer than $maximum_length bp, and all sequences that were present less than $minimum_count times were discarded." >> ${dms}_results/read_me.txt
fi 
