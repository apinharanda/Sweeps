#!/bin/bash

#this script takes the single fasta done from the msms_to_fasta
#adds uncertainty of haps for outgroup letters
#modifies the names of the samples
#pastes the sequences to be of the desired lenght
#run in the same dir as the fastas that have been converted from msms to fasta

#v1 9June 2020
#v2 Oct 2022

#run this script after the check_sweeps_fixed.sh (splits msms output and checks which fixed)
#and after msms_to_fasta.py (makes the msms output into letter fastas and incorporates uncertainty as random numbers)
#it works to paste any number of sequences together, and should be used even when not pasting anything together

#0) help -h or --help . If yes display usage
if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then
        echo "Usage: bash $0 [list of all msms.fasta to paste / format for ASF]"
        exit 0
fi


#usage
#bash make_multi_fasta_from_msms_fasta.sh

args=("$@")

echo "PART 3a    Make multi fasta & format for ASF"
echo "bash multi_fasta_msms_v2.sh $@"


echo Number of msms-fastas to paste: $#

echo "List of msms-fastas:"
for var in "$@"
do
    echo "$var"
done

rm -rf *tmp

###read from command line
#create tmp fastas

for var in "$@"; do
	cat "$var" > $var".tmp" 
done

#paste sequences together
echo "List of msms-fastas tmp:"

for var in *.tmp; do
	echo "$var"
done

#make name out fasta
#concatenate args
echo "Making msms multi fasta......"

function concatenate_args
{
    string=""
    for a in "$@" # Loop over arguments
    do
        if [[ "${a:0:1}" != "-" ]] # Ignore flags (first character is -)
        then
            if [[ "$string" != "" ]]
            then
                string+="_" # Delimeter
            fi
            string+="$a"
        fi
    done
    echo "$string"
}

#paste_give name of concatenated string
#in the order user gave it
#format so that it can go directly to AFS
#change numbers in outgroup with uncertainty to letters for normal fasta
#paste sequences in order, remove spaces, remove extra names >

paste "$@".tmp | sed 's/ ,//g' | sed 's/\[//g' | sed 's/\]//g' | sed 's/,//g' | perl -pe 's/\S+/abs($&) > 0.0225 ? "A" : $&/ge' | sed -E 's/0.[0-9]+/T/g' |   sed -e "s/[[:space:]]\+//g" | sed 's/>/\t>/g' | awk '{print $1}' | sed 's/[0-9]*//g' >  "$(concatenate_args "$@")"


echo "Name of msms multi fasta: $(concatenate_args "$@")"


#####Make input for ELS
####
#3b) then have the counts for each class 
echo "PART 3b    Count derived mutations and format input for ELS"
echo "perl count_derived_mutations_ELS_v2.pl "$(concatenate_args "$@") > "$(concatenate_args "$@").count.log" 

perl count_derived_mutations_ELS_v2.pl $(concatenate_args "$@") > $(concatenate_args "$@").count.log

mv ELS_results $(concatenate_args "$@").ELS_results


###then this second part creates individual fastas for pi analysis
#this creates separate fastas by greping the name +A1


#remove tmp files
rm -rf *.tmp


