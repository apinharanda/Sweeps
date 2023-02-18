#-t: 4*1,000,000*100,000*2.8e-9
#-r: 4*1,000,000*100,000*2.46e-8
#-I: npop, size1, size2, mrate
#-Sc: t, deme, 2Ns for AA, 2Ns for Aa, 2Ns for aa
#-SI: t (in 4Ne), ndeme, initial allele freq of beneficial allele in deme1: 1/2/1,000,000, in deme2
#-ej: 2,000,000(generations)/4/1,000,000

nohup java -Xmx128G -jar ../msms/lib/msms.jar 18 1 -t 11200 -r 98400 -seed 0x59594c2b8e99cbdc -I 2 10 8 0 -N 1000000 -Sc 0 1 200000 100000 0 -SI 0.4999998 2 5e-07 0 -ej 0.5 2 1 -Sp 0.5 -Smark > test3.txt &

#check genotype at selected site
#awk 'NR==6{for(i=1;i<=NF;i++){if($i=="0.50000"){midind=i;};};}NR>=7{split($0,hap,"");print hap[midind-1];}' test3.txt
#awk 'hapcount>0{split($0,hap,"");fixedder*=hap[midind];hapcount--;if (hapcount==0 && fixedder==1){print simnum;}}NR%22==6{fixedder=1;simnum++;for(i=1;i<=NF;i++){if($i=="0.50000"){midind=i;};};hapcount=8;}' test2.txt
#awk 'NR==23,NR==46' test1.txt > test2.txt

python msms_to_fasta.py

#check pi in each deme
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < test2.fasta

for file in *.fa 
do
	sed -i '1c>test2'  ${file}
done

module load gcc/4.9.3 
./calculatePolymorphism Dsan*.fa | ./nonOverlappingWindows -o Dsan_pi_w5000.tsv -w 5000 -n 
./calculatePolymorphism Dyak*.fa | ./nonOverlappingWindows -o Dyak_pi_w5000.tsv -w 5000 -n
./calculatePolymorphism D*.fa | ./nonOverlappingWindows -o Dsan_Dyak_pi_w5000.tsv -w 5000 -n 

python SDM_v4.py
