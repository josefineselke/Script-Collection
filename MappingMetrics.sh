#!/usr/bin/env bash

# this script was used in my master thesis to summarize
# metrics from read mappings of simulated datasets

tools=("Art" "neat" "simuscop")
samples=("v42_axicel" "v42_R11")

for sample in ${samples[@]}; do
  for tool in ${tools[@]}; do
    files=${sample}/${tool}*/*filtered.sam
    all=$(grep -c ^@ /path/to/directory/${tool}*/${sample}/*.fq)
    calcReads=$(wc -l </path/to/directory/${tool}*/${sample}/*.fastq)

    for file in ${files[@]}; do

      mapped=$(wc -l <$file)
      countTP=0
      countFN=0
      filtered=/path/to/directory/${tool}*/${sample}/*.fastq

      while read line; do
        TP=$(grep -c "${line:1}"$'\t' $file)
        if [ "$TP" = 0 ]; then
          let countFN++
        else
          let countTP=countTP+$TP
        fi
      done < $filtered

      let FP=$mapped-$countTP

      AllMappedSAM="${file%.*}"
      AllMappedFiltered="${AllMappedSAM%.*}.sam"
      AllMapped=$(grep -v '^@' ${AllMappedFiltered} \
                | awk '!($2 == 4)' | wc -l)
      echo -e "${file}\t${AllMapped}" >> AllMapped.txt
      let TN=${Allmapped}-${countTP}-${countFN}-${FP}

      echo -e "${file}\t${Allmapped}\t${mapped}\t${countTP} \
             \t${FP}\t${countFN}\t${TN}\t${all}\t${calcReads}" \
      >> metrics.txt
    
    done
  done
done
