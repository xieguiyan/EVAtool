#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at dot hust dot edu dot cn
# @DATE: 2021-09-03 11:46:34
# @DESCRIPTION:

from pathlib import Path
from config import Config
from logger import Logger
import subprocess


class Fastq(object):
    def __init__(self, inputfile: Path, outputdir: Path, config: Config, logger: Logger):
        self.inputfile = Path(inputfile)
        self.outputdir = outputdir
        self.config = config
        self.logger = logger
        self.trimname = f"{self.inputfile.stem}.fastq.filter.1.gz"
        # self.processfq = self.process_fastq()

    def is_sra(self):
        return True if self.inputfile.suffix == ".sra" else False

    def dump_fastq(self):
        cmd = [self.config.config["fastqdump"], "--split-files", self.inputfile]
        print(cmd)
        return subprocess.run(cmd)

    def trim(self) -> None:
        sra_filter_reads_1 = self.trimname
        filter_params = f"ILLUMINACLIP:{self.config.config['adp_path']}:{self.config.config['trimmomatic_sRNA_para']}"
        cmd = f"java -jar -Xms8000m -Xmx8000m {self.config.config['trimmomatic']} SE -threads {self.config.config['cpu_number']} {self.inputfile.stem}_1.fastq.gz {sra_filter_reads_1} {filter_params} -trimlog {sra_filter_reads_1}.log 1>{sra_filter_reads_1}_run.log 2>&1"
        return subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")

    def process_fastq(self):
        if self.is_sra():
            rc = self.dump_fastq()
            if rc.returncode == 0:
                self.logger(f"Dump SRA file {self.inputfile} to fastq")
            else:
                self.logger.logtofile(f"Error in dump SRA file {self.inputfile} to fastq")

        runtrim = self.trim()
        if runtrim == 0:
            self.logger.logtofile(f"Trimm {self.inputfile.stem} fastq file")
        else:
            self.logger.logtofile(f"Error in trimm {self.inputfile.stem} fastq file")


fastq = Fastq(
    inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra",
    outputdir="/home/xiegy/github/EVAtool/test/example-data",
    config=Config(),
    logger=Logger(message="preocessed fastq!", logfile="/home/xiegy/github/EVAtool/src/evatool/utils/log.txt"),
)

# fastq.process_fastq()
# print(fastq.trimname)
# print(fastq.outputdir)
# print(fastq.config.config)
print(fastq.logger.logtofile)
# fastq.outputdir / fastq.trimname
