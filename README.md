# **EVAtool**

---
###### Background
non-coding RNAs (ncRNA) in extracellular vesicles (EVs) role as important agents of cell-to-cell communication. The abundance of transcripts varies in different states, especially in the tumors. Therefore, quantification of ncRNAs in EVs plays an important role in the process of exploring cancer diagnostic biomarkers.

###### ++Requires++
python >= 3.5 <br>

###### Softwares
samtools = 1.12 <br>
bowtie2 = 2.4.2 <br>
fastq-dump = 2.10.9 <br>
trimmomatic-0.39.jar <br>
bedtools = 2.30.0 <br>

###### Usage
- Reference and config file download
    - wget "http://bioinfo.life.hust.edu.cn/EVAtool/ref/refs.zip"
    - The details of Reference:
        - Homo_sapiens.GRCh38.86.chr.gtf
        - Homo_sapiens_v86
        - miRNA miRBase V22
        - snoRNA snoDB V1.2.1
        - tRNA GtRNAdb V18.1
        - piRNA RNAcentral V17
        - rRNA RNAcentral V17
        - snRNA NCBI GenBank
        - YRNA NCBI GenBank

    - The content of config file:
        -   "trimmomatic_sRNA_para": "2:10:4:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15",
        -   "tag_cut": "0",
        -   "cpu_number": "8",
        -   "RPM": "total",
        -   "bowtie_para_4_miRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 2 -L 15 -i S,1,1.15",
        -   "bowtie_para_4_rRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15",
        -   "bowtie_para_4_tRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15",
        -   "bowtie_para_4_piRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15",
        -   "bowtie_para_4_snRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 3 -L 10 -i S,1,1.15",
        -   "bowtie_para_4_snoRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -L 12 -i S,1,1.15",
        -   "bowtie_para_4_YRNA": "--no-head --no-sq -D 15 -R 2 -N 1 -k 1 -L 12 -i S,1,1.15",
        -   "bowtie_para_4_genome": "-D 15 -R 2 -N 1 -k 1 -L 10 -i S,1,1.15",
        -   "adp_path": "/refs/sRNA.fa",
        -   "miRNA_index": "/refs/miRNA/hsa.hairpin.fa",
        -   "mirbase": "/refs/miRNA/hsa.mirbase.txt",
        -   "piRNA_index": "/refs/piRNA/piRNA.fa",
        -   "YRNA_index": "/refs/YRNA/YRNA.fa",
        -   "snoRNA_index": "/refs/snoRNA/snoRNA.fa",
        -   "snRNA_index": "/refs/snRNA/snRNA.fa",
        -   "rRNA_index": "/refs/rRNA/rRNA.fa",
        -   "tRNA_index": "/refs/tRNA/tRNA.fa",
        -   "genome_index": "/refs/Homo_sapiens_v86/Homo_sapiens_v86.chromosomes.fasta",
        -   "gtf_annotation": "/refs/Homo_sapiens.GRCh38.86.chr.gtf",
        -   "trans_bed_annotation": "/refs/four_elements.sort.bed",
        -   "exon_bed_annotation": "/refs/two_elements.sort.bed"

- Install evatool
   
```
pip install -i https://test.pypi.org/simple/ evatoolt
```

- Example

```
>>>from evatool.utils.fastq import Fastq
>>>from evatool.utils.config import Config
>>>from evatool.utils.logger import Logger
>>> def test_config():
...     from evatool.utils.config import Config
...     config = Config(configfile="./refs/reference_config.json")
...     print(config.config)

>>> test_config()
```
- Input parameters

```
#-i: sra file with path (required)
#-o: output directory (required)
#-c: configure file with path ( not required)
#-n: ncRNA type list (not required)
```

- Two ways of testing the scripts

```
bash /home/xiegy/github/EVAtool/test/example-script/example_evatool.sh

or

/home/xiegy/github/EVAtool/venv/bin/python \
/home/xiegy/github/EVAtool/evatool/main.py \
-i /home/xiegy/github/EVAtool/test/example-data/example.fastq.gz \
-o /home/xiegy/github/EVAtool/test/tmp_result/example_fq_gz \
#-c /home/xiegy/github/EVAtool/evatool/resource/configure.json \
#-n "miRNA" "rRNA" "tRNA" "piRNA" "snoRNA" "snRNA" "YRNA"
```

###### Modules in evatool
1. Logger
    - Produce a log file that records the events when the evatool is running.
    - Logger(logfile = Path)
2. Config
    - Tansfer the content of config file to each module. 
    - Config(config = Path)
3. Fastq
    - Transform the raw .sra file to .fastq or .fastq.gz, quality control (QC) and adapter trimming.
    -  Fastq(inputfile, outputdir, config=Config, log=Logger, ncrna_lst=ncrna_lst)
4. Tag
    - Statistics the distribution of various types of reads and transformed the trimmed file into .fasta.
    - Tag(fastq=Fastq(inputfile, outputdir, config=Config, log=Logger, ncrna_lst=ncrna_lst))
5. SAM
    - Map reads to human genome and other types ncRNA reference.
    - Tag(fastq=Fastq(inputfile, outputdir, config=Config, log=Logger, ncrna_lst=ncrna_lst))
6. Bam
    - Transform sam file to bam file.
    - Bam(fastq=Fastq, tag=Tag)
7. Stat
    - Quantification and RPM normalization.
    - Stat(fastq=Fastq, tag=Tag)
8. Plot
    - Visualize output results.
    - Plot(inputfile, outputdir)
9. Report
    - Provide online results report.
    - Report(inputfile, outputdir, config=Config, ncrna_lst=ncrna_lst))


###### Advanced usage
1. Use each module independently
>     After installing evatool, users only need to import the interested package(s) in the python interpreter to use the corresponding module.
> e.g. 

```
>>>from evatool.utils.fastq import Fastq
>>>from evatool.utils.config import Config
>>>from evatool.utils.logger import Logger
>>> def test_config():
...     from evatool.utils.config import Config
...     config = Config(configfile="./refs/reference_config.json")
...     print(config.config)

>>> test_config()
```

2. Custom ncRNA type(s) 
    - Based on existing reference sequences
        - change the input parameter of ncRNA type list, the default ncRNAs are include 7 types: "miRNA" "rRNA" "tRNA" "piRNA" "snoRNA" "snRNA" "YRNA".
    - Add other type reference sequences
        - Three steps :
        
        1. Add the ncRNA reference and index into refs;

        ```
        mkdir [ncRNA name]
        add the reference and index to the directory
        ```

        2. Change the input parameter of ncRNA type list;
        
        ```
        -n  "miRNA" "rRNA" "tRNA"
        ```

        3. Add the ncRNA name and reference directory in config file as following the existing ways.
        
        ```
        "miRNA_index": "/refs/miRNA/hsa.hairpin.fa"
        ```



###### Directory tree (main)
├── __init__.py <br>
├── main.py <br>
├── resource <br>
│   ├── __init__.py <br>
│   ├── logging.conf <br>
│   ├── template_report.html <br>
│   └── tool_config.json <br>
└── utils <br>
    ├── bam.py <br>
    ├── config.py <br>
    ├── fastq.py <br>
    ├── __init__.py <br>
    ├── logger.py <br>
    ├── plot.py <br>
    ├── README.md <br>
    ├── report.py <br>
    ├── sam.py <br>
    ├── stat.py <br>
    └── tag.py <br>


###### Last news
-  evatoolt 0.1.1 upload evatool to pypi
-  evatoolt 0.1.2 add README.md
 
 