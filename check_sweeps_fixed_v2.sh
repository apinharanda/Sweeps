#!/bin/bash

#find which simulations from msms output fixed
#v1 26May2020
#v2 Oct12022


#run where the log or out file with all the collated simunations from msms is 
#usage

#0) help -h or --help . If yes display usage
if [[ ( $@ == "--help") ||  $@ == "-h" ]]
then 
	echo "Usage: $0 [name dir for output] [out file from msms run] [position of selection] [number indv simulated per pop]"
	exit 0
fi 

#bash check_sweeps_fixed.sh recent_5000_0.006 recent_5000_0.006.log 0.50000000 20
#$1 name of folder recent_5000_0.006
#$2 name of outfile recent_5000_0.006.log
#$3 site where selection should occur 0.50000000
#$4 how many individuals of each

##### read in command line
folder=$1
name_of_outfile=$2
site=$3
number_samples_per_pop=$4

echo "PART 1    Checking which sweeps fixed"
echo "For neutral simulations the dir $1_simulations_fixed should be empty"
echo "bash check_sweeps_fixed.sh $1 $2 $3 $4"

#make dir like the start of the output msms file
mkdir $1

cd $1

#remove double // from the end of each simuation
sed -i 's/\/\///g' ../$2

# split the large file with all the simulations into individual ones each time it encounters segsites
#give the same name as the main file 
csplit -f $1"_" -z ../$2 /segsites/ '{*}'

# remove the 00 as its empty
rm -rf $1"_00"

# find the column of genotypes for the position that we decided it was under selection
for i in $1"_"*; do
    awk -v var=$3 'NR==2{for(i=1;i<=NF;i++){if($i==var){midind=i;};};}NR>=3{split($0,hap,"");print hap[midind-1];}' $i | sed -e "s/[[:space:]]\+/\t/g" > $i.tmp
done

#every number of rows as we have individuals in the population sum
#same sums + file name in temp file 2

for i in $1"_"*".tmp"; do
	awk -v var=$4 '{s+=$1}NR%var==0{print s;t+=s;s=0}END{print "total_fixed: ",t,FILENAME;if(s) print "left: " s}' $i > $i.tmp2
done

#create folder
mkdir ../$1"_"simulations_fixed

#those temp2 that have a sum equal or higher than the number of individuals for either the first n or the second n rows
#this means the sweep fixed
#save in file 

for i in  $1"_"*".tmp2"; do
	cat $i | tr "\n" "\t" | awk -v var=$4 '{if(($1 >= var) || ($2 >= var)) print $5}' | sed 's/.tmp//g'  >> simulations_fixed.log
done

#open file and move the simulations to the new folder
for file in $(cat simulations_fixed.log); do mv "$file" ../$1"_"simulations_fixed; done

#create folder for the logs of the fixed
mkdir ../$1"_"simulations_fixed_logs

#then open the file again and move the tmps too
for file in $(cat simulations_fixed.log); do mv "$file.tmp" ../$1"_"simulations_fixed_logs; done
for file in $(cat simulations_fixed.log); do mv "$file.tmp.tmp2" ../$1"_"simulations_fixed_logs; done

#remove all the rest that didnt fix
#rm -rf ../$1

#remove tmp files
rm *.tmp*

#move fixed logs to fixed log dir
mv ../$1/simulations_fixed.log ../$1"_"simulations_fixed_logs

#change the dir name of the ones that didnt fix 
mv ../$1 ../$1"_notFixed"


