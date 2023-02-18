#!/bin/bash

#AP bash usage of msms_fasta
#fasta out var

#0) help -h or --help . If yes display usage
if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then 
	echo "Usage: bash $0 [length of simulated sequence] [total nr of haps simulated] [output msms after check_sweeps_fixed.sh, full path to indv file] [speciesID 1] [speciesID 2] [nr haps species 1] [nr haps species 2]"
	exit 0
fi 

#size_sequence = $1 #1000000
#number_haps = $2 #40
#msms_folder_name_file = $3 #the invididual fixed or not file after running check_sweeps_fixed.sh
#sp_1 = $4 #Dsim
#sp_2 = $5 #Dmau
#outgroup = $6 #outgroup_Dmel
#number_haps_1 = $7 #20
#number_haps_2 = $8 #20

echo "PART 2    Changing msms output into fasta format"
echo "python msms_to_fasta_v3.py $1 $2 $3 $4 $5 $6 ${7} ${8} $3.fasta" 

python msms_to_fasta_v3.py $1 $2 $3 $4 $5 $6 ${7} ${8} $3.fasta

