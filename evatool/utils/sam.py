#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-26 11:06:09
# @DESCRIPTION:

from .fastq import Fastq
from .config import Config
import subprocess
from pathlib import Path


class SAM(object):
    def __init__(self, fastq: Fastq):
        self.fastq = fastq
        self.samdir = Path(fastq.inputfile)
        self.ncrna_lst = ["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "scRNA"]

    def ncRNA_map(self) -> None:
        sample_input = self.fastq.inputfile.stem
        tag_fa = f"{self.fastq.outputdir}/{self.fastq.inputfile.stem}.fa"
        for i in self.ncrna_lst:
            output_sam = f"{self.fastq.outputdir}/{sample_input}.{i}.sam"
            bowtie_stat = f"{self.fastq.outputdir}/{sample_input}.{i}.bowtie.stat"
            cmd = [
                self.fastq.config.config["bowtie"],
                "-x",
                self.fastq.config.config[i + "_index"],
                "-p",
                self.fastq.config.config["cpu_number"],
                self.fastq.config.config["bowtie_para_4_{0}".format(i)],
                "-f",
                tag_fa,
                "-S",
                output_sam,
                "2>",
                bowtie_stat,
            ]
            cmd = " ".join(cmd)
            print(cmd)
            map_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if map_result.returncode == 0:
                self.fastq.log.log(f"Sucess in align to {i} reference!")
            else:
                self.fastq.log.log(f"Failed in align to {i} reference!")
