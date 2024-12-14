#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at dot hust dot edu dot cn
# @DATE: 2021-09-03 11:46:34
# @DESCRIPTION:

import os
import sys

from pathlib import Path
import subprocess

from .config import Config
from .logger import Logger

class Fastq(object):
    def __init__(self, inputfile: Path, outputdir: Path, config: Config, log: Logger, ncrna_lst: list = ["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "YRNA"]):
        self.inputfile = Path(inputfile)
        self.outputdir = outputdir
        self.config = config
        self.log = log
        self.ncrna_lst = ncrna_lst
        self.trimname = f"{self.inputfile.stem}.fastq.trimmed.gz"
        
        

    def is_sra(self):
        return True if self.inputfile.suffix == ".sra" else False

    def is_fastq(self):
        return True if (self.inputfile.suffix == ".fastq" or self.inputfile.suffix == ".fq" or "".join(self.inputfile.suffixes) == ".fastq.gz" or "".join(self.inputfile.suffixes) == ".fq.gz") else False

    def dump_fastq(self):
        cmd = [self.config.config["fastqdump"], "--dumpbase", "--gzip", "--split-files", "-O", self.outputdir, self.inputfile]
        return subprocess.run(cmd)
    
    def detect_quality_format(self,file_path):
        command = (
            f"zcat {file_path} | head -100 | awk '{{if(NR%4==0) printf(\"%s\",$0);}}' | "
            "od -A n -t u1 | awk 'BEGIN{min=500;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{"
            "if(max<=73 && min<59) print \"-phred33\"; "
            "else if(max<=104 && min>=64) print \"-phred64\"; "
            "else if(min>=59 && min<64 && max>73) print \"Solexa+64\"; "
            "else print \"-phred33\";}'"
        )
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
        return result.stdout.strip()
    
    def trim(self, fq_file) -> None:
        trimmed_fq = f"{self.outputdir}/{self.trimname}"
        filter_params = f"ILLUMINACLIP:{self.config.config['adp_path']}:{self.config.config['trimmomatic_sRNA_para']}"
        cmd = f"{self.config.config['trimmomatic']} SE {self.quality_format} -threads {self.config.config['cpu_number']} {fq_file} {trimmed_fq} {filter_params}"        
        print(cmd)
        return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
    
    def process_fastq(self):
        if self.is_sra():
            rc = self.dump_fastq()
            if rc.returncode == 0:
                fq_file = str(f"{self.outputdir}/{self.inputfile.stem}_1.fastq.gz")
                self.quality_format = self.detect_quality_format(fq_file)
                runtrim = self.trim(fq_file)
                if runtrim == 0:
                    self.log.log(message=f"Success in trimm {self.inputfile.stem} fastq file")
                else:
                    self.log.log(message=f"Error in trimm {self.inputfile.stem} fastq file")
                    sys.exit(1)
                self.log.log(message=f"Success in dump SRA file {self.inputfile.stem} to fastq")
            else:
                self.log.log(message=f"Error in dump SRA file {self.inputfile.stem} to fastq")
                sys.exit(1)
        elif self.is_fastq():
            fq_file = self.inputfile
            self.quality_format = self.detect_quality_format(fq_file)
            runtrim = self.trim(fq_file)
            if runtrim.returncode == 0:
                self.log.log(message=f"Success in trimm {self.inputfile.stem} fastq file")
            else:
                self.log.log(message=f"Error in trimm {self.inputfile.stem} fastq file")
                sys.exit(1)
        