#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Gui
# @CONTACT: xieguiyan.at.hust.dot.edu.dot.com
# @DATE: 2021-09-03 11:21:07
# @DESCRIPTION:

import sys
from pathlib import Path

sys.path.append("../EVAtool")

from evatool.utils.fastq import Fastq
from evatool.utils.config import Config
from evatool.utils.logger import Logger
from evatool.utils.sam import SAM
from evatool.utils.bam import Bam
from evatool.utils.tag import Tag
from evatool.utils.stat import Stat
from evatool.utils.report import Report

import argparse


def run(inputfile: Path, outputdir: Path, config: Path, ncrna_lst: list) -> None:
    config = Config(config)
    logger = Logger(f"{outputdir}.log.txt")
    fastq_result = Fastq(inputfile, outputdir, config=config, log=logger, ncrna_lst=ncrna_lst)
    fastq_result.process_fastq()
    tag_result = Tag(fastq=fastq_result)
    tag_result.pocess_stat()
    sam_result = SAM(fastq=fastq_result)
    sam_result.get_sam()
    bam_result = Bam(fastq=fastq_result, tag=tag_result)
    bam_result.process_bam()
    stat_result = Stat(fastq=fastq_result, tag=tag_result)
    stat_result.stat_match()
    return 1


def main(configure):
    print(configure)
    result = run(configure.input, configure.output, configure.config, configure.ncrna)
    if result == 1:
        report_result = Report(configure.input, configure.output)
        report_result.prepare_html()
        print("Success!")
    else:
        print("Failed!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="EVAtool could be used to estimate the quantification and abundence of small ncRNA from EV or other sources.",
    )
    parser.add_argument("-i", "--input", help="The path of the input file, the file type could be '.sra, .fastq.gz or .fastq'.", required=True)
    parser.add_argument("-o", "--output", help="The path of output file.", required=True)
    parser.add_argument("-c", "--config", help="The path of the Config file. User can download the config file from url, or define yourself.", required=False)
    parser.add_argument("-n", "--ncrna", nargs="*", type=str, help="The list of small ncRNA types.  User can use default ncRNA list (miRNA, rRNA, tRNA, piRNA, snoRNA, snRNA, YRNA), or define yourself.", required=False)
    configure = parser.parse_args()
    main(configure)
