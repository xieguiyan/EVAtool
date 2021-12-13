# EVAtool

---
###### Background
non-coding RNAs (ncRNA) in extracellular vesicles (EVs) role as important agents of cell-to-cell communication. The abundance of transcripts varies in different states, especially in the tumors. Therefore, quantification of ncRNAs in EVs plays an important role in the process of exploring cancer diagnostic biomarkers.
###### Requires
python >= 3.5 <br>
###### Softwares
samtools = 1.12 <br>
bowtie = 2.4.2 <br>
fastq-dump = 2.10.9 <br>
trimmomatic-0.39.jar <br>
bedtools = 2.30.0 <br>

###### Usage
- Reference and config file download
    - wget "xxxxxx"

- install evatool
   
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
- Command line parameters

```
#-i: sra file with path (required)
#-o: output directory (required)
#-c: configure file with path ( not required)
#-n: ncRNA type list (not required)
```

- Two ways of running the scripts

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

###### Configure
1. xxx
2. xxx
3. xxx


###### Advanced usage
1. xxx
2. xxx
3. xxx

###### Directory tree
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


###### 版本内容更新
 evatoolt 0.1.1 add README.md
 
 