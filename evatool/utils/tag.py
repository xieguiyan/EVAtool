#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-07 10:19:59
# @DESCRIPTION:

from pathlib import Path
from .fastq import Fastq
import gzip
import subprocess

class Tag(object):
    def __init__(self, fastq: Fastq):
        self.fastq = fastq
        self.prefix = f"{self.fastq.outputdir}/{self.fastq.inputfile.stem}"
        self.freqfile = f"{self.prefix}.freq.stat"
        self.tagfile = f"{self.prefix}.fa"
        self.tag_count_dict = self.get_tag_count()

    def is_trimm(self):
        trimmfile = Path(f"{self.fastq.outputdir}/{self.fastq.trimname}")
        return trimmfile.exists()

    def stat_tag(self):
        tag_dict = {}
        f = gzip.open(self.fastq_output_path)
        for n, line in enumerate(f.readlines()):
            line_number = n + 1
            if line_number % 4 != 2:
                continue
            fq_seq = line.strip()
            try:
                tag_dict[fq_seq] += 1
            except KeyError:
                tag_dict[fq_seq] = 1
        sorted_tag_number = sorted(tag_dict.keys(), key=lambda z: tag_dict[z], reverse=True)
        return sorted_tag_number, tag_dict

    def store_tag(self, sorted_tag_number, tag_dict):
        fq_len_frequency_dict = {}
        tag_count = {}
        reads_n = 0
        out_reads_n = 0
        with open(self.tagfile, "w") as tf:
            for n, line in enumerate(sorted_tag_number):
                seq_len = len(line)
                seq_count = tag_dict[line]
                if seq_len in fq_len_frequency_dict:
                    fq_len_frequency_dict[seq_len] += seq_count
                else:
                    fq_len_frequency_dict[seq_len] = seq_count
                reads_n += seq_count
                tag_number = n + 1
                tag_number = "t{0:0>8d}".format(tag_number)
                cut_off = int(self.fastq.config.config["tag_cut"])
                if seq_count > cut_off:
                    out_reads_n += seq_count
                    tf.write(f">{tag_number}\t{seq_count:d}\n{line}\n")
                    # tag_count[tag_number] = int(seq_count)
        return fq_len_frequency_dict, reads_n, out_reads_n

    def store_freq(self, fq_len_frequency_dict, reads_n, out_reads_n):
        len_stat_lst = [0, 0, 0, 0]
        with open(self.freqfile, "w") as of:
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
                of.write(f"{i}\t{len_n}\t{len_freq:f}\n")
            flag = "ok" if len_stat_lst[1] + len_stat_lst[2] > 0.5 else "no"
            of.write("{0}\t{1}\t{2:d}\t{3:.2f}\t{4:d}\n".format(flag, "\t".join([str("{0:.2f}".format(j)) for j in len_stat_lst]), out_reads_n, out_reads_n / reads_n, reads_n))

    def get_tag_count(self):
        if Path(self.tagfile).exists():
            tag_count = {}
            with open(self.tagfile, "r") as f:
                for i in f:
                    if i.startswith(">"):
                        tag_info = i.strip(">\n").split("\t")
                        tag_count[tag_info[0]] = int(tag_info[1])
            return tag_count

    def pocess_stat(self):
        if self.is_trimm():
            self.fastq_input_path = Path(self.fastq.outputdir) / self.fastq.trimname
            self.fastq_output_path = Path(self.fastq.outputdir) / f"qc_{self.fastq.trimname}"
            quality_format = "" if self.fastq.quality_format == "-phred33" else "--phred64"
            fastp_command = [
                "fastp",
                quality_format,
                "--n_base_limit", "5",
                "--length_required", "15",
                "--length_limit", "50",
                "--average_qual", "20",
                "--disable_adapter_trimming",
                "--thread", self.fastq.config.config['cpu_number'],
                "--json", f"{self.prefix}.qc_report.json",
                "--html", f"{self.prefix}.qc_report.html",
                "-i", str(self.fastq_input_path),
                "-o", str(self.fastq_output_path),
            ]
            try:
                # 执行fastp命令
                print(fastp_command)
                subprocess.run(fastp_command, check=True)
                self.fastq.log.log(message=f"Fastp QC completed: {self.fastq_output_path}")
            except subprocess.CalledProcessError as e:
                self.fastq.log.log(message=f"Fastp QC failed: {e}")
                return 
            
            sorted_tag_number, tag_dict = self.stat_tag()
            fq_len_frequency_dict, reads_n, out_reads_n = self.store_tag(sorted_tag_number, tag_dict)
            self.store_freq(fq_len_frequency_dict, reads_n, out_reads_n)
            if Path(self.tagfile).exists():
                self.tag_count_dict = self.get_tag_count()
                self.fastq.log.log(message="Success in stat seq in fq!")
            else:
                self.fastq.log.log(message="Error in stat seq in fq!")
        else:
            self.fastq.log.log(message="Trimmed file is not exist!")
