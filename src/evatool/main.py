#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Gui
# @CONTACT: xieguiyan.at.hust.dot.edu.dot.com
# @DATE: 2021-09-03 11:21:07
# @DESCRIPTION:

import os, sys
from typing import Any
from pathlib import Path
from evatool.utils import Fastq
from evatool.utils import Config
from evatool.utils import Logger
import argparse


def run(inputfile: Path, config: Config, outputdir: Path) -> None:
    config = Config(config).config
    logger = Logger(f"{outputdir}.log.txt")
    fastqfile = Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sh", outputdir="/home/xiegy/github/EVAtool/test/example-data")
    print(config)
    print("Hello World!")


def main():
    run()


if __name__ == "__main__":
    main()
