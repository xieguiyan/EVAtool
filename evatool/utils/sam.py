#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-26 11:06:09
# @DESCRIPTION:

from evatool.utils.fastq import Fastq
from .stat import Stat
from .logger import Logger
from .config import Config
import subprocess
from pathlib import Path


class SAM(object):
    def __init__(self, fastq: Fastq) -> None:
        self.fastq = fastq
        self.samdir = Path(fastq.inputfile)

    def ncRNA_map(self, ncRNA_lst) -> None:
        sample_input = self.fastq.inputfile.stem
        tag_fa = f"{self.fastq.outputdir}/{self.fastq.inputfile.stem}.fa"
        for n, i in enumerate(ncRNA_lst):
            locals()["{0}_sam".format(i)] = "{0}.{1}.sam".format(sample_input, i)
            locals()["{0}_bowtie_stat".format(i)] = "{0}.{1}.bowtie.stat".format(sample_input, i)
            cmd = [
                self.config.config["bowtie"],
                "-x",
                self.config.config[i + "_index"],
                "-p",
                self.config.config["cpu_number"],
                self.config.config["bowtie_para_4_{0}".format(i)],
                "-f",
                tag_fa,
                "-S",
                locals()["{1}/{0}_sam".format(i, Path(self.fastq.outputdir))],
                "2>",
                locals()["{1}/{0}_bowtie_stat".format(i, Path(self.fastq.outputdir))],
            ]
            cmd = " ".join(cmd)
            map_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if map_result.returncode == 0:
                self.logger.log(f"Sucess in align to {i} reference!")
            else:
                self.logger.log(f"Failed in align to {i} reference!")
