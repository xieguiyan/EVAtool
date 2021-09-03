#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 11:46:34
# @DESCRIPTION:

import os
from pathlib import Path
from evatool.utils.config import Config
import sys
import datetime
import logging
import json
import subprocess

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a file handler
handler = logging.FileHandler("benchmarks.log")
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
# add the file handler to the logger
logger.addHandler(handler)


class Fastq(object):
    """Fastq class"""

    def __init__(self, inputfile: Path, config: Config, outputdir: Path):
        self.inputfile = inputfile
        self.workdir = outputdir
        self.config = config

    def __str__(self):
        return "SRA file: {}\nWork diretory: {}\nConfig: {}".format(self.inputfile, self.workdir, self.config_file)

    def is_sra(self):
        return True if self.inputfile.suffix == ".sra" else False

    def dump_fastq(self):
        cmd = [self.config["fastqdump"], "--split-files", self.inputfile]
        return subprocess.check_output(args=cmd)

    def trim(self) -> None:
        pass

    def process_fastq(self):
        if self.is_sra():
            rc = self.dump_fastq()
            if rc == 0:
                logger.info("Dump SRA file {} to fastq".format(self.inputfile))
            else:
                logger.info("Error in dump SRA file {} to fastq".format(self.inputfile))

        self.trim()

    def transfq(self):
        resultdir = "{dir}/analysis.bowtie.{time}".format(dir=self.workdir, time=datetime.datetime.now().strftime("%Y-%m-%d"))
        isExist = os.path.isdir(resultdir)
        if not isExist:
            os.makedirs(resultdir)
        os.chdir(resultdir)
        with open(self.config_file, "r") as c:
            config_env = json.load(c)
        logger.info("Directory has been created!\n")
        sra = os.path.splitext(os.path.basename(self.inputfile))[0]
        sra_reads = "{sra}_1.fastq.gz".format(sra=sra)
        quals = "--dumpbase"
        transfq = "{fastqdump} {quals} --gzip --split-files {inputfile} 1>{sra}.dump.log 2>&1".format(fastqdump=config_env["fastqdump"], quals=quals, inputfile=self.inputfile, sra=sra)
        runfq = subprocess.run(transfq, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        if runfq.returncode == 0:
            logger.info("sucess in transfer sra to fastq!")
        else:
            logger.info("error in transfer sra to fastq!")
        sra_filter_reads_1 = "{sra}.fastq.filter.1.gz".format(sra=sra)
        filter_params = "ILLUMINACLIP:{adpter}:{filter_params}".format(adpter=config_env["adp_path"], filter_params=config_env["trimmomatic_sRNA_para"])
        trimfq = "java -jar -Xms8000m -Xmx8000m {trimmomatic} SE -threads {cpu_number} {sra_reads} {out_reads} {filter_params} -trimlog {out_reads}.log 1>{out_reads}_run.log 2>&1".format(
            trimmomatic=config_env["trimmomatic"],
            cpu_number=config_env["cpu_number"],
            sra_reads=sra_reads,
            out_reads=sra_filter_reads_1,
            filter_params=filter_params,
        )
        runtrim = subprocess.run(trimfq, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        if runtrim.returncode == 0:
            logger.info("sucess in trimmomatic for fastq!")
        else:
            logger.info("error in trimmomatic for fastq!")
        return "{}\n".format(resultdir)
