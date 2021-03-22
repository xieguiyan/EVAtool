cd /home/xiegy/github/EVAtool/test/example-data/DRR006758
mkdir -p /home/xiegy/github/EVAtool/test/example-data/DRR006758/analysis.bowtie.4pq8fqia.2021-03-23 && cd /home/xiegy/github/EVAtool/test/example-data/DRR006758/analysis.bowtie.4pq8fqia.2021-03-23
echo 'start processing' > benchmarks.log && date >> benchmarks.log
echo 'start fastq-dump' >> benchmarks.log && date >> benchmarks.log
/home/xiegy/github/EVANEQ/bin/fastq-dump  --gzip --split-files /home/xiegy/github/EVAtool/test/example-data/SRR10078125.sra 1>SRR10078125.dump.log 2>&1
java -jar -Xms8000m -Xmx8000m /home/xiegy/github/EVANEQ/bin/trimmomatic-0.36.jar SE -threads 8 SRR10078125_1.fastq.gz SRR10078125.fastq.filter.1.gz ILLUMINACLIP:/home/xiegy/github/EVANEQ/refs/sRNA.fa:2:10:4:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 -trimlog SRR10078125.fastq.filter.1.gz.log 1>SRR10078125.fastq.filter.1.gz_run.log 2>&1
java -jar -Xms8000m -Xmx8000m /home/xiegy/github/EVANEQ/bin/trimmomatic-0.36.jar SE -threads 8 SRR10078125.fastq.filter.1.gz SRR10078125.fastq.filter.fianl.gz ILLUMINACLIP:/home/xiegy/github/EVANEQ/refs/sRNA.fa:2:10:4:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 -trimlog SRR10078125.fastq.filter.fianl.gz.log 1>SRR10078125.fastq.filter.fianl.gz_run.log 2>&1
echo 'end fastq-dump' >> benchmarks.log && date >> benchmarks.log
echo 'start reads2tag translation' >> benchmarks.log && date >> benchmarks.log
