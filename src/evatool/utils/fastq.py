#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at dot hust dot edu dot cn
# @DATE: 2021-09-03 11:46:34
# @DESCRIPTION:

from pathlib import Path
from evatool.utils.config import Config
from evatool.utils.logger import Logger
import subprocess


class Fastq(object):
    """Fastq class"""

    def __init__(self, inputfile: Path, config: Config, outputdir: Path):
        self.inputfile = inputfile
        self.outputdir = outputdir
        self.config = config
        self.logger = Logger()

    def is_sra(self):
        return True if self.inputfile.suffix == ".sra" else False

    def dump_fastq(self):
        cmd = [self.config["fastqdump"], "--dumpbase --split-files", self.inputfile]
        return subprocess.check_output(args=cmd)

    def trim(self) -> None:
        sra_filter_reads_1 = f"{self.inputfile.stem}.fastq.filter.1.gz"
        filter_params = f"ILLUMINACLIP:{self.config['adp_path']}:{self.config['trimmomatic_sRNA_para']}"
        cmd = f"java -jar -Xms8000m -Xmx8000m {self.config['trimmomatic']} SE -threads {self.config['cpu_number']} {self.inputfile.stem}_1.fastq.gz {sra_filter_reads_1} {filter_params} -trimlog {sra_filter_reads_1}.log 1>{sra_filter_reads_1}_run.log 2>&1"
        runtrim = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")

    def process_fastq(self):
        if self.is_sra():
            rc = self.dump_fastq()
            if rc == 0:
                Logger.mylogger(f"Dump SRA file {self.inputfile} to fastq")
            else:
                Logger.mylogger(f"Error in dump SRA file {self.inputfile} to fastq")

        runtrim = self.trim()
        if runtrim == 0:
            Logger.mylogger(f"Trimm {self.inputfile.stem} fastq file")
        else:
            Logger.mylogger(f"Error in trimm {self.inputfile.stem} fastq file")
