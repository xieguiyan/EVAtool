#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-07 10:19:59
# @DESCRIPTION:

from pathlib import Path
import pathlib
from fastq import Fastq


class Stat(object):
    def __init__(self, statfile: Path, fastq: Fastq):
        self.statfile = statfile
        self.fastq = fastq

    def is_sra(self):
        trimmfile = pathlib.Path(self.fastq.trimname)
        return trimmfile.exists()

    def stat_tag(self):
        tag_dict = {}
        file_type = "txt"
        len_stat_lst = [0, 0, 0, 0]
        i = self.fastq.trimname
        file_t = commands.getoutput("file %s" % (i))
        if re.search(r"gzip compressed data", file_t):
            file_type = "gz"
        elif re.search(r"ASCII text", file_t):
            file_type = "txt"
        if file_type == "gz":
            f = gzip.open(i, "rb")
        elif file_type == "txt":
            f = open(i, "r")
        for n, i in enumerate(f.readlines()):
            line_number = n + 1
            if line_number % 4 != 2:
                continue
            fq_seq = i.strip()
            try:
                tag_dict[fq_seq] += 1
            except KeyError:
                tag_dict[fq_seq] = 1
        sorted_tag_number = sorted(tag_dict.keys(), key=lambda z: tag_dict[z], reverse=True)
        return sorted_tag_number

    def store_tag(self):
        fq_len_frequency_dict = {}
        reads_n = 0
        out_reads_n = 0
        with open(f"{self.fastq.inputfile.stem}.fa", "w") as of:
            for n, i in enumerate(sorted_tag_number):
                fq_len = len(i)
                tag_count = tag_dict[i]
                if fq_len in fq_len_frequency_dict:
                    fq_len_frequency_dict[fq_len] += tag_count
                else:
                    fq_len_frequency_dict[fq_len] = tag_count
                reads_n += tag_count
                tag_number = n + 1
                tag_number = "t{0:0>8d}".format(tag_number)
                if tag_count > self.fastq.config["tag_cut"]:
                    out_reads_n += tag_count
                    of.write(">{0}\t{1:d}\n{2}\n".format(tag_number, tag_count, i))

    def stat_freq(self):
        with open(configure["outfile"] + ".freq.stat", "w") as of:
            for i in fq_len_frequency_dict:
                len_n = fq_len_frequency_dict[i]
                len_freq = len_n / reads_n
                if i < 15:
                    len_stat_lst[0] += len_freq
                elif 15 <= i < 30:
                    len_stat_lst[1] += len_freq
                elif 30 <= i <= 40:
                    len_stat_lst[2] += len_freq
                else:
                    len_stat_lst[3] += len_freq
                of.write("{0}\t{1}\t{2:f}\n".format(i, len_n, len_freq))
            flag = "ok" if len_stat_lst[1] + len_stat_lst[2] > 0.5 else "no"
            of.write("{0}\t{1}\t{2:d}\t{3:.2f}\t{4:d}\n".format(flag, "\t".join([str("{0:.2f}".format(j)) for j in len_stat_lst]), out_reads_n, out_reads_n / reads_n, reads_n))
