#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-26 11:06:09
# @DESCRIPTION:

from sys import argv
from .fastq import Fastq
from .config import Config
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count


class SAM(object):
    def __init__(self, fastq: Fastq):
        self.fastq = fastq
        self.samdir = Path(fastq.inputfile)

    def ncRNA_map(self, ncrna) -> None:
        sample_input = self.fastq.inputfile.stem
        tag_fa = f"{self.fastq.outputdir}/{self.fastq.inputfile.stem}.fa"
        output_sam = f"{self.fastq.outputdir}/{sample_input}.{ncrna}.sam"
        bowtie_stat = f"{self.fastq.outputdir}/{sample_input}.{ncrna}.bowtie.stat"
        cmd = [
            self.fastq.config.config["bowtie"],
            "-x",
            self.fastq.config.config[ncrna + "_index"],
            "-p",
            self.fastq.config.config["cpu_number"],
            self.fastq.config.config["bowtie_para_4_{0}".format(ncrna)],
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
            self.fastq.log.log(f"Sucess in align to {ncrna} reference!")
        else:
            self.fastq.log.log(f"Failed in align to {ncrna} reference!")

    def get_sam(self) -> None:
        ncrna_lst = ["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "scRNA"]
        p = Pool(cpu_count())
        for ncrna in ncrna_lst:
            print(ncrna)
            p.apply_async(self.ncRNA_map, args=(ncrna,))
        p.close()
        p.join()
