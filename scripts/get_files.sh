#!/bin/bash
Final_proj="/athena/angsd/scratch/mac4017/Final_Project"
cd $Final_proj
samples_to_download=( $(egrep "RNASeq" $Final_proj/metadata/elife-30860-supp1-v1.tsv | cut -d$'\t' -f 3) )
samp_names=( $(egrep "RNASeq" $Final_proj/metadata/elife-30860-supp1-v1.tsv | cut -d$'\t' -f 1) )
mkdir Fastqs
cd Fastqs
mkdir ${samp_names[@]}

for i in ${!samples_to_download[@]}; do
   cd ${samp_names[$i]}
   echo downloading ${samp_names[$i]} fastqs: $i of ${#samples_to_download[@]}
   fastqs=$(egrep  ${samples_to_download[$i]} $Final_proj/metadata/filereport_read_run_PRJEB12982_tsv.txt | cut -d$'\t' -f9)
   IFS=' ' read -r -a fqs <<< $fastqs
   for x in ${fqs[@]}; do
      IFS=";" read -r -a p_fq <<< $x
      wget -t 2 -nc -w 30 ftp://${p_fq[0]}
      sleep 10
      wget -t 2 -nc -w 30 ftp://${p_fq[1]}
      sleep 15
   done
   for x in *; do
      mv $x ${samp_names[$i]}_${x}
   done
   cd ..
done
