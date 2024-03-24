#!/usr/bin/env bash

# this script was used in my master thesis to use various quantification software for scRNA-seq data more easily

# kallisto | bustools: kb-python/0.27.3-foss-2021b , can process gzip
# Salmon/1.9.0-GCC-11.2.0, can process gzip
# STAR/2.7.9a-GCC-11.2.0, can process gzip
# CellRanger/6.1.2, can process gzip

while getopts 't:s:p:r:w:o:c:n:h:a:m' OPTION; do
  case "$OPTION" in
    a)
      threads=$OPTARG
    m)
      mem=$OPTARG
    t)
      tool=$OPTARG
      ;;
    s)
      sample=$OPTARG
      ;;
    p)
      path=$OPTARG
      ;;
    r)
      reference=$OPTARG
      ;;
    w)
      whitelist=$OPTARG
      ;;
    o)
      output=$OPTARG
      ;;
    c)
      chemistry=$OPTARG
      ;;
    n)
      num_cells=$OPTARG
      num_cells="--expect-cells ${num_cells}"
      ;;
    h)
      echo "Usage:"
      echo "-a: Number of threads available for analysis"
      echo "-m: Available RAM in GB"
      echo "-t: Tool to Use [alevin, kallisto, STARsolo, cellranger]"
      echo "-s: Samples to Analyze"
      echo "-p: Path to Data in CellRanger-readable format, without / at the end"
      echo "-r: Reference to align to; all necessary files need to be in an index/ directory in the project"
      echo "-w: Whitelist"
      echo "-o: Output directory"
      echo "-c: chemistry (10XV2, 10XV33, 10XV25, SC5P-R2)"
      echo "-n: Number of expected cells"
      echo "Provide all arguments in quotation marks separated by spaces!"
      ;;
    ?)
      echo "Use -h to get help."
      exit 0
      ;;
  esac
done
shift "$(($OPTIND -1))"

for t in ${tool}; do
  tools+=("$t")
done
for s in ${sample}; do
  samples+=("$s")
done
for r in ${reference}; do
  references+=("$r")
done
for p in ${path}; do
  path+=("$p")
done

if [ "${chemistry}" == "10XV2" ]; then
  a_tech=chromium
  kb_tech=10XV2
  S_UMI=10
elif [ "${chemistry}" == "10XV33" ]; then
  a_tech=chromiumV3
  kb_tech=10XV3
  S_UMI=12
elif [ "${chemistry}" == "10XV25" ]; then
  a_tech=chromium
  kb_tech=10XV2
  S_UMI=10
elif [ "${chemistry}" == "SC5P-R2" ]; then
  cr_tech="--chemistry SC5P-R2"
  a_tech=chromium
  kb_tech=10XV2
  S_UMI=10

fi

if [ ! -d "/path/to/project/${output}/" ]; then
  mkdir /path/to/project/${output}
fi

for sample in ${samples[@]}; do
  cd /path/to/project
  if [ ! -d "/path/to/project/${output}/${sample}/" ]; then
    mkdir /path/to/project/${output}/${sample}
  fi
  for tool in ${tools[@]}; do
    if [ ! -d "/path/to/project/${output}/${sample}/${tool}/" ]; then
      mkdir /path/to/project/${output}/${sample}/${tool}
    fi
  done


  readUMI=$(ls -1 ${path[1]}/${sample}*L001_R1*)
  read=$(ls -1 ${path[1]}/${sample}*L001_R2*)
  readUMI2=$(ls -1 ${path[1]}/${sample}*L002_R1*)
  read2=$(ls -1 ${path[1]}/${sample}*L002_R2*)


    # Salmon-Alevin
    if [[ " ${tools[*]} " =~ " alevin " ]]; then
      if [ ! -d "/path/to/project/${output}/${sample}/alevin/${ref}" ]; then
        mkdir /path/to/project/${output}/${sample}/alevin/${ref}
      fi

      salmon alevin \
      -l ISR \
      -1 ${readUMI} ${readUMI2}\
      -2 ${read} ${read2}\
      --${a_tech} \
      -i index/alevin/${ref}/ \
      -p ${threads} \
      -o path/to/project/${output}/${sample}/alevin/${ref} \
      --tgMap index/${ref}.tsv
    fi

    # Kallisto-Bustools
    if [[ " ${tools[*]} " =~ " kallisto " ]]; then
      if [ ! -d "path/to/project/${output}/${sample}/kallisto/${ref}" ]; then
        mkdir path/to/project/${output}/${sample}/kallisto/${ref}
      fi

      kb count \
      -t ${threads} \
      --h5ad \
      -x ${kb_tech} \
      -i index/kallisto/${ref}/${ref}.idx \
      -g index/kallisto/${ref}/${ref}.txt \
      -o path/to/project/${output}/${sample}/kallisto/${ref}/ \
      --filter bustools \
      ${readUMI} \
      ${read} \
      ${readUMI2} \
      ${read2}


    fi

    # STARsolo
    if [[ " ${tools[*]} " =~ " STARsolo " ]]; then
      if [ ! -d "path/to/project/${output}/${sample}/STARsolo/${ref}" ]; then
        mkdir path/to/project/${output}/${sample}/STARsolo/${ref}
      fi

      STAR \
      --runThreadN ${threads} \
      --genomeDir index/star/${ref}/ \
      --soloType Droplet \
      --readFilesCommand gunzip -c \
      --outFilterMultimapNmax 1 \
      --soloUMIlen ${S_UMI} \
      --soloCBwhitelist $whitelist \
      --outFileNamePrefix quantification/${output}/${sample}/STARsolo/${ref}/quant.${ref}. \
      --readFilesIn ${read},${read2} ${readUMI},${readUMI2}
    fi

    # CellRanger
    if [[ " ${tools[*]} " =~ " cellranger " ]]; then
      if [ ! -d "path/to/project/${output}/${sample}/cellranger/${ref}" ]; then
        mkdir path/to/project/${output}/${sample}/cellranger/${ref}
      fi

      cd /path/to/project/${output}/${sample}/cellranger/${ref}
      cellranger count \
      --id=${sample} ${cr_tech} \
      --fastqs=${path}/ \
      --transcriptome=/path/to/project/index/cellranger/${ref}/ \
      --sample=${sample} \
      --nosecondary \
      --localcores=${threads} \
      --localmem=${mem} ${num_cells}

    fi
  done
done
