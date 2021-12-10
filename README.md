# EVAtool
Extracellular Vesicles Abundance of ncRNAs expression quantification tool

[Home page](http://bioinfo.life.hust.edu.cn/EVAtool/#)

#bin files for softwares
1. samtools = 1.12
2. bowtie = 2.4.2
3. fastq-dump = 2.11.0
5. trimmomatic-0.39.jar
6. bedtools = 2.30.0

Usage:

#-i: sra file with path (required)
#-o: output directory (required)
#-c: configure file with path ( not required)
#-n: ncRNA type list (not required)

#bash /home/xiegy/github/EVAtool/test/example-script/example_evatool.sh

# or

/home/xiegy/github/EVAtool/venv/bin/python \
  /home/xiegy/github/EVAtool/evatool/main.py \
  -i /home/xiegy/github/EVAtool/test/example-data/example.fastq.gz \
  -o /home/xiegy/github/EVAtool/test/tmp_result/example_fq_gz
  #-c /home/xiegy/github/EVAtool/evatool/resource/configure.json \
  #-n "miRNA" "rRNA" "tRNA" "piRNA" "snoRNA" "snRNA" "YRNA"
  

![work-flow](https://user-images.githubusercontent.com/19505178/130918235-e4fde4bc-5e6b-4450-a7e1-39f2427256eb.png)

