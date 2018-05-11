#! /bin/sh
# -*- coding: utf-8 -*-
#SBATCH --cpus-per-task 4
#SBATCH -p batch
#SBATCH --mem 20G

# First use mkproject 
# then load the module for the project
# then run this script

# Path variables
project="${project_root}"
short_reads="${raw}"
project_code="${project_name}"
log="${logs}/qc.log"

# Program parameters
cpus=4
minlen=52
qscore=5

ffqs=($(find ${short_reads} -type l -regex ".*\.forward\..*\.gz" | sort))
rfqs=($(find ${short_reads} -type l -regex ".*\.reverse\..*\.gz" | sort))
for ((i=0; i<${#ffqs[*]}; i++)); do
  forward=${ffqs[i]};
  reverse=${rfqs[i]};
  
  # Left and right reads are pairs
  IFS=". " read -a right <<< $(basename ${forward});
  IFS=". " read -a left <<< $(basename ${reverse});
  sample=${right[0]};
  if [[ "$sample" != "${left[0]}" ]]; then
     echo "$(date +"%b %-d %k:%M:%S") QC failed. Read pairs are out of order." >> ${log};
     break;
  fi
 
  # Raw read statistics
  rawstats="${stats}/${project_code}.raw_stats.csv";
  if [[ ! -f $rawstats ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: generating basic statistics on raw reads." >> ${log};
    readstats.py --csv --output ${rawstats} ${forward} ${reverse};
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: finished generating basic statistics on raw reads." >> ${log};
  else
    rawinfo=$(grep -i "$sample" $rawstats);
    if [ -z "$rawinfo" ]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: generating basic statistics on raw reads." >> ${log};
      readstats.py --csv ${forward} ${reverse} | tail -n +2 >> ${rawstats};
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: finished generating basic statistics on raw reads." >> ${log};
    else
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: raw read statistics have already been generated ... skipping." >> ${log};
    fi
  fi

  # Read QC
  qc_paired="${qc}/${sample}.interleaved.atrim.decontam.qtrim.fq.gz";
  qc_singles="${qc}/${sample}.singles.atrim.decontam.qtrim.fq.gz";
  if [[ ! -f $qc_paired ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: starting quality control of sequencing reads." >> ${log};
    bbduk.sh -Xmx6g threads=2 qin=33 overwrite=t ref=/srv/databases/contaminants/truseq_adapters.fa in1=${forward} in2=${reverse} out=stdout.fq stats=${stats}/${sample}.adapt_stats.txt ftm=5 ktrim=r k=23 mink=9 rcomp=t hdist=2 tbo tpe minlength=0 2>${logs}/${sample}.adapters.log | bbduk.sh -Xmx8g threads=1 qin=33 overwrite=t interleaved=t in=stdin.fq out=stdout.fq outm=${qc}/discarded/${sample}.phix.fq.gz ref=/srv/databases/contaminants/phix174.fa.gz k=31 hdist=1 mcf=0.9 stats=${stats}/${sample}.phix_stats.txt 2>${logs}/${sample}.decontam.log | bbduk.sh -Xmx2g ziplevel=9 threads=1 qin=33 overwrite=t interleaved=t in=stdin.fq out=${qc_paired} outs=${qc_singles} stats=${stats}/${sample}.trim_stats.txt qtrim=rl trimq=${qscore} minlength=${minlen} 2>${logs}/${sample}.qtrim.log;
    sleep 30s;
    if [[ -s $qc_paired ]]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: quality control of sequencing reads completed." >> ${log};
    else
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: quality control of sequencing reads failed. ${qc_paired} is empty." >> ${log};
      continue
    fi
  else
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: quality control of sequencing reads has already been performed ... skipping." >> ${log};
  fi

  # Final read statistics
  qcstats="${stats}/${project_code}.qc_stats.csv";
  if [[ ! -f $qcstats ]]; then
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: generating basic statistics on quality-controlled short reads." >> ${log};
    readstats.py --csv --output ${qcstats} ${qc_paired};
    echo "$(date +"%b %-d %k:%M:%S") ${sample}: finished generating basic statistics on quality-controlled short reads." >> ${log};
  else
    qcinfo=$(grep -i "$sample" $qcstats);
    if [ -z "$qcinfo" ]; then
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: generating basic statistics on quality-controlled short reads." >> ${log};
      readstats.py --csv ${qc_paired} | tail -n +2 >>${qcstats}
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: finished generating basic statistics on quality-controlled short reads." >> ${log};
    else
      echo "$(date +"%b %-d %k:%M:%S") ${sample}: basic read statistics have already been generated on quality-controlled short reads ... skipping." >> ${log};
    fi
  fi
done
