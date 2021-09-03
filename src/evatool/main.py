#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 11:21:07
# @DESCRIPTION:

import os, sys
from src.evatool.utils import fastq
from pathlib import Path
from evatool.utils import Fastq
from evatool.utils import Config
from evatool.utils import Logger


def run(dir: Path, inputfile, outputdir, configfile) -> None:
    config = Config(configfile)
    logger = Logger(dir / f"{outputdir}.log.txt")
    fastqfile = Fastq("a", "b", "c")
    fastqfile.process_fastq()
    a.test()
    print("Hello World!")


def main():
    run()


if __name__ == "__main__":
    main()
