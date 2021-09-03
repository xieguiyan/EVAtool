#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 11:21:07
# @DESCRIPTION:

import os, sys
from evatool.utils import Fastq


def run():
    a = Fastq("a", "b", "c")
    a.test()
    print("Hello World!")


def main():
    run()


if __name__ == "__main__":
    main()
