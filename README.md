# Sweeps
Scripts to run several sweep finding pipelines - this is very much work in progress


msms_to_fasta.py script originally wrote by Clair Han; count_derived_mutations_ELS_v2.pl by Peter Andolfatto

Follow ms-ASF.log on how to go from msms output to fasta file

Part 1, msms format output is ready to be used in a streamlined way.

The steps are:

——msms + format output

1) run msms

2) check_sweeps_fixed.sh 
check in which simulations there was fixation and where there wasn’t 
this is important for weak/recent selection as most don't fix
It serves as a sanity check for the rest, and it does process the msms output (splits it into separate runs)

3) msms_to_fasta_v1.sh
split msms output to fasta file
run in python2 and with module np installed

4) multi_fasta_msms_v2.sh 
combines x number of simulations into one concatenated multi-fasta file formatted for ASF pipeline
counts derived mutations for ELS input 


Part 2, ASF run still needs some work, but the ground work is done.

——ASF

5) calculate pi

6) genome_derived_expectation.py

7) normalization.R

8) SDM_v7.py

9) Plot

