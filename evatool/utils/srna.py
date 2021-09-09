#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-09-08 20:10:26
# @DESCRIPTION:

from .stat import Stat
from .logger import Logger
from .config import Config
import re
import subprocess


class SRNA(object):
    def __init__(self, stat: Stat, logger: Logger, config: Config):
        self.stat = stat
        self.logger = logger
        self.config = config

    def deal_mir_info(self):
        mature_miRNA = {}
        hairpin_info = {}
        with open(self.config.config["mirbase"], "r") as f:
            for i in f:
                line = i.strip().split("\t")
                hairpin_seq = line[4]
                mir_seq = line[3]
                miRNA_start = hairpin_seq.find(mir_seq)
                if miRNA_start == -1:
                    continue
                miRNA_len = len(mir_seq)
                hairpin_len = len(hairpin_seq)
                arm = "5p" if hairpin_len > 2 * miRNA_start + miRNA_len - 1 else "3p"
                mature_miRNA.setdefault(line[1], {})[arm] = [line[0], miRNA_start, miRNA_len, mir_seq, hairpin_seq]
                if line[1] in hairpin_info:
                    continue
                hairpin_info[line[1]] = [hairpin_seq, hairpin_len]
        return (mature_miRNA, hairpin_info)

    def get_fa4sRNA(self):
        tag_count = {}
        with open(self.stat.tagfile, "r") as f:
            for line in f:
                if line.startswith(">"):
                    tag_info = line.strip(">\n").split("\t")
                    tag_count[tag_info[0]] = int(tag_info[1])
        return tag_count

    def get_true_miRexp(self, tag_hairpin_dict, hairpin_tag_dict, mature_miRNA, tag_count_dict):
        divide_tag_expressoin_dict = {}
        mir_exp_dict = {}
        tag_with_multi_assign_dict = {tag: ref for tag, ref in tag_hairpin_dict.items() if len(ref) > 1}
        sorted_tagn_4_multi_maps = sorted(tag_with_multi_assign_dict.keys(), key=lambda a: tag_count_dict[a], reverse=True)
        for tag in tag_hairpin_dict:
            hairpins = list(tag_hairpin_dict[tag].keys())
            hairpins_n = len(hairpins)
            tag_count = tag_count_dict[tag]
            refs_mir = set()
            tmp_dict = {}
            for hairpin in hairpins:
                mapped_details = tag_hairpin_dict[tag][hairpin]
                for mapped_detail in mapped_details:
                    arm = mapped_detail[-1]
                    if arm in mature_miRNA[hairpin]:
                        mir = mature_miRNA[hairpin][arm][0]
                    else:
                        arm = str(8 - int(arm[0])) + "p"
                        mir = hairpin + "-" + arm
                    refs_mir.add(mir)
                    tmp_dict.setdefault(mir, set()).add(hairpin)
            mir_n = len(refs_mir)
            if mir_n == 1:
                tag_count_division = tag_count
            else:
                tag_count_division = int(tag_count / mir_n) + 1
            for j in refs_mir:
                if j in mir_exp_dict:
                    mir_exp_dict[j] += tag_count_division
                else:
                    mir_exp_dict[j] = tag_count
        return mir_exp_dict

    def ncRNA_map(self, ncRNA_lst) -> None:
        sample_input = self.stat.fastq.inputfile
        tag_fa = self.stat.tagfile
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
                locals()["{0}_sam".format(i)],
                "2>",
                locals()["{0}_bowtie_stat".format(i)],
            ]
            cmd = " ".join(cmd)
            map_result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            if map_result.returncode == 0:
                self.logger.log(f"Sucess in align to {i} database!")
            else:
                self.logger.log(f"Failed in align to {i} database!")

    def ncRNA_exp(self, ncRNA_lst):
        mapped_nc_tag_dict = {}
        ref_exp = {}
        miR_mapped_detail = {}
        mapped_ncRNA_counts = {}
        distance_re = re.compile(r"NM:i:(\d)")
        mismatch_pattern = re.compile(r"MD:Z:(\d+)[A-Z]?(\d*)")
        ref_tag_detail = {}
        tag_ref_detail = {}
        (mature_miRNA, hairpin_info) = self.deal_mir_info(self.config.config["mirbase"])
        for n, i in enumerate(ncRNA_lst):
            ref_tag_detail[i] = {}
            tag_ref_detail[i] = {}
            with open(locals()["{0}_sam".format(i)], "r") as sf:
                for j in sf:
                    line = j.strip().split("\t")
                    if line[1] != "0":
                        continue
                    (tag_name, mapping_detail, mapped_ref) = line[0:3]
                    m = distance_re.search(line[-1])
                    if m:
                        mapped_distance = int(m.group(1))
                    else:
                        mapped_distance = 0
                    if tag_name.startswith("@") or mapped_distance > 1:
                        continue
                    if n:
                        if tag_name in mapped_nc_tag_dict:
                            continue
                        tag_ref_detail[i].setdefault(tag_name, set()).add(mapped_ref)
                    else:
                        arm = "3p"
                        if hairpin_info[line[2]][-1] - int(line[3]) * 2 - len(line[9]) + 1 > 0:
                            arm = "5p"
                        mapped_start = int(line[3])
                        if mapped_distance:
                            discard_flag = 0
                            mismatch_detail = mismatch_pattern.search(line[-2])
                            mismatch_loci = int(mismatch_detail.group(1)) + mapped_start - 1
                            #                        mature_miRNA.setdefault(line[1],{})[arm] = [line[0],mapped_start,miRNA_len]
                            mir_seed_regions = [(z[1] + 1, z[1] + 7) for z in mature_miRNA[line[2]].values()]
                            for region in mir_seed_regions:
                                if region[0] <= mismatch_loci <= region[1]:
                                    discard_flag = 1
                            if discard_flag:
                                continue
                        try:
                            tag_ref_detail[i].setdefault(tag_name, {})[mapped_ref].append([mapped_start, mapped_distance, arm])
                        except KeyError:
                            tag_ref_detail[i].setdefault(tag_name, {})[mapped_ref] = [[mapped_start, mapped_distance, arm]]
                    ref_tag_detail[i].setdefault(mapped_ref, set()).add(tag_name)
                    mapped_nc_tag_dict[tag_name] = i

        tag_count_dict = self.get_fa4sRNA(tag_fa)
        for n, i in enumerate(ncRNA_lst):
            ref_exp[i] = {}
            mapped_tags = tag_ref_detail[i].keys()
            mapped_ncRNA_counts[i] = sum([tag_count_dict[tag] for tag in mapped_tags])
            if n:
                for ref in ref_tag_detail[i]:
                    ref_counts = sum([tag_count_dict[xx] for xx in ref_tag_detail[i][ref]])
                    ref_exp[i][ref] = ref_counts
            else:
                tags_mir_count = {}
                arm_lst = set(["5p", "3p"])
                mir_exp_dict = self.get_true_miRexp(tag_ref_detail[i], ref_tag_detail[i], mature_miRNA, tag_count_dict)
                ref_exp[i] = mir_exp_dict

    return (ref_exp, mapped_ncRNA_counts, tag_count_dict, mapped_nc_tag_dict, ref_tag_detail)

    def process_map(self):
        ncRNA_lst = ["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "scRNA"]
        (ref_exp, mapped_ncRNA_counts, tag_count_dict, mapped_nc_tag_dict, ref_tag_detail) = self.ncRNA_exp(ncRNA_lst)
