#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 11:46:34
# @DESCRIPTION:


class Fastq(object):
    """Fastq class"""

    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "@{}\n{}\n+\n{}\n".format(self.name, self.seq, self.qual)

    def test(self):
        print(self.__str__())
