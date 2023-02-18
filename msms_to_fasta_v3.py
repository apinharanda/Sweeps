#clairh script
#AP added sys vars 2Jun20 incorporated misreference if the ancestral state into the simulations
#Oct2022 added usage

#usage
#python msms_to_fasta.py seq_length num_haplo filename species1 species2 outgroup hap1 hap2 output_file
#python msms_to_fasta.py 5000 40 recent_5000_0.006.log Dsim Dmau outgroup_Dmel 20 20 recent_5000_0.006.fasta

#seq_length lenght of simulated sequence
#num_haplo how many total haps were simulated
#filename needs to be the parsed output from msms, one file per sequence, use check_sweeps_fixed.sh first to only use the ones that fixed for selection
#species1
#species2
#outgroup
#hap1 number of haps of species 1
#hap2 number of haps of species 2
#output_file name of output fasta
 
 
import numpy 
import sys
import random
import numpy as np
import argparse

seq_length=int(sys.argv[1])
num_haplo=int(sys.argv[2])

filename=sys.argv[3]
specie1=sys.argv[4]
specie2=sys.argv[5]
outgroup=sys.argv[6]
hap1=int(sys.argv[7])
hap2=int(sys.argv[8])
output_file=sys.argv[9]


filereader=open(filename)
output_fasta=open(output_file,"w")

segsites=None
positions=[]
genotypes=None

counter=0
for line in filereader:
	if counter<num_haplo:
		if segsites!=None and len(positions)>0:
			genotypes[counter,:]=numpy.array(list(line.rstrip()))
			counter+=1
		if line[0:8]=="segsites":
			segsites=int(line.rstrip().split(' ')[1])
			genotypes=numpy.zeros(shape=(num_haplo,segsites), dtype=bool)
		if line[0:9]=="positions":
			positions=numpy.round(numpy.array(line.rstrip().split(' ')[1:],dtype=float)*seq_length)

genotypes_full = numpy.zeros((num_haplo,seq_length), dtype=bool)  
for i in range(0, segsites):
	genotypes_full[:,int(positions[i])-1]=genotypes[:,i]

output = numpy.chararray((num_haplo,seq_length)) 
output[:] = "A" 
output[genotypes_full] = "T"

#round(random.uniform(greaterThan, lessThan), digits)
ran_floats = [round(random.uniform(0,1),4) for i in xrange(1000000)]

for i in range(0, hap1):
	output_fasta.write('>' + specie1 + str(i+1) + '\n'+ ''.join(output[i,:]) + '\n')
for i in range(0, hap2):
	output_fasta.write('>' + specie2 + str(i+1) + '\n'+ ''.join(output[i+hap1,:]) + '\n')
###need to incorporate misreference of the ancestral state
output_fasta.write('>' + outgroup + '\n'+ ''.join([str(ran_floats)]) + '\n')


output_fasta.close()


