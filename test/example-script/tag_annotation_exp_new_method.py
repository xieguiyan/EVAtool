#! /usr/bin/env python
"""
@@ script: sra_proceduere_bowtie.py
@@ description: based on python3.5.2
@@ author: zhangq
@@ modified: xiegy
"""


from collections import Counter
from operator import itemgetter
import argparse
import datetime
import os
import sys
import re
import subprocess
import tempfile
import copy
import itertools
import pandas


def printlg(log_message):
    print(log_message, file=log_file)


def get_fa4sRNA(tag_fa):
    tag_count = {}
    with open(tag_fa, "r") as f:
        for i in f:
            if i.startswith(">"):
                tag_info = i.strip(">\n").split("\t")
                tag_count[tag_info[0]] = int(tag_info[1])
    return tag_count


def get_config_env(config_file):
    config_env = {}
    with open(config_file) as cf:
        for i in cf:
            if i.startswith("#"):
                continue
            if i.find("=") > 1:
                fields = [j.strip() for j in i.strip().split("=")]
                config_env[fields[0]] = fields[1]
    return config_env


def deal_mir_info(mirbase_file):
    mature_miRNA = {}
    hairpin_info = {}
    with open(mirbase_file, "r") as f:
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


def get_true_miRexp(tag_hairpin_dict, hairpin_tag_dict, mature_miRNA, tag_count_dict):
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


def get_tag_flag(stat_file):
    cmd = ["tail", "-n", "1", stat_file]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = res.stdout.decode().split("\n")
    fields = lines[0].strip().split("\t")
    if fields[0] == "ok":
        kept_tags = int(fields[-3])
        return kept_tags
    else:
        printlg("The length of tags were out of the range 15-40 nt, please check the data. Time: {0}".format(datetime.datetime.now().strftime("%c")))
        sys.exit()


def sub_dict_modify(dct, keys, flag):
    if flag == 1:
        return dict([(key, dct.get(key)) for key in keys])
    else:
        return dict([(key, dct.get(key)) for key in dct.keys() - keys])


def deal_mapped_info(mapped_anno_bed, mapped_unanno_bed, mapped_nc_tag_dict, un_annotated_tags_nc, ref_exp, tag_count_dict):
    details_mapped_tag = {}
    map_to_genome_tags = set()
    new_anno_tag_detail = {}
    region_anno_detail = {}
    rebuilt_unanno_info = []
    tmp_anno = []
    un_annotated_tags = copy.deepcopy(un_annotated_tags_nc)
    with open(mapped_unanno_bed, "r") as uf:
        for i in uf:
            line = i.strip().split("\t")
            tags_id = line[4].strip().split(",")
            map_to_genome_tags.update(tags_id)
            tags_in_nc = [tag_id for tag_id in tags_id if tag_id in mapped_nc_tag_dict]
            tag_type_lst = [mapped_nc_tag_dict[tag] for tag in tags_in_nc]
            if tags_in_nc:
                tag_type_counter = Counter(tag_type_lst)
                sorted_tag_type_counter = sorted(tag_type_counter.keys(), key=lambda a: tag_type_counter[a], reverse=True)
                tag_type = sorted_tag_type_counter[0]
                if "other" in ref_exp[tag_type]:
                    ref_exp[tag_type]["other"] += sum([tag_count_dict[t_z] for t_z in tags_id if t_z not in mapped_nc_tag_dict])
                else:
                    ref_exp[tag_type]["other"] = sum([tag_count_dict[t_z] for t_z in tags_id if t_z not in mapped_nc_tag_dict])
                tmp_anno.extend(tags_id)
                for t_z in tags_id:
                    mapped_nc_tag_dict[t_z] = tag_type
            else:
                rebuilt_unanno_info.append(i)
        un_annotated_tags = sub_dict_modify(un_annotated_tags, set(tmp_anno), 2)
    with open(mapped_anno_bed, "r") as mf:
        anno_c_set = []
        new_anno_tags = set()
        un_annotated_tags_tmp = un_annotated_tags.keys()
        for i in mf:
            line = i.strip().split("\t")
            ref_gene_info = line[-9:]
            mapped_info = line[0:7]
            tags_id = mapped_info[4].strip().split(",")
            map_to_genome_tags.update(tags_id)
            unanno_flag = [tag_id for tag_id in tags_id if tag_id in un_annotated_tags_tmp]
            tags_in_nc = [tag_id for tag_id in tags_id if tag_id in mapped_nc_tag_dict]
            if unanno_flag:
                if tags_in_nc and tags_in_nc not in anno_c_set:
                    tag_type_lst = [mapped_nc_tag_dict[tag] for tag in tags_in_nc]
                    tag_type_counter = Counter(tag_type_lst)
                    sorted_tag_type_counter = sorted(tag_type_counter.keys(), key=lambda a: tag_type_counter[a], reverse=True)
                    tag_type = sorted_tag_type_counter[0]
                    if "other" in ref_exp[tag_type]:
                        ref_exp[tag_type]["other"] += sum([tag_count_dict[t_z] for t_z in tags_id if t_z not in mapped_nc_tag_dict])
                    else:
                        ref_exp[tag_type]["other"] = sum([tag_count_dict[t_z] for t_z in tags_id if t_z not in mapped_nc_tag_dict])
                    new_anno_tags.update(unanno_flag)
                    for t_z in tags_id:
                        mapped_nc_tag_dict[t_z] = tag_type
                    anno_c_set.append(tags_in_nc)
            strand = "genebody" if mapped_info[3] == ref_gene_info[5] else "opposite"
            try:
                details_mapped_tag[mapped_info[4]][strand].append(ref_gene_info)
            except KeyError:
                details_mapped_tag.setdefault(mapped_info[4], {})[strand] = [ref_gene_info]
                details_mapped_tag[mapped_info[4]]["l"] = mapped_info
        #                un_annotated_tags = sub_dict_modify(un_annotated_tags, tags_id, 2)
        un_annotated_tags = sub_dict_modify(un_annotated_tags, new_anno_tags, 2)
    for tags in details_mapped_tag:
        strand_lst = ["genebody", "opposite"]
        tag_anno_dict = dict([(z, "-") for z in strand_lst])
        tag_ids = tags.strip().split(",")
        # tag_mapped_catagory_handle.write('TagId\tTagCount\tCategory\tmapped_item\tregion\n')
        detail_info = details_mapped_tag[tags]
        mapped_info = detail_info["l"]
        mapped_start = int(mapped_info[1])
        mapped_end = int(mapped_info[2])
        for strand_flag in strand_lst:
            if strand_flag not in detail_info:
                continue
            ref_lst = detail_info[strand_flag]
            try:
                tmp_lst = sorted(ref_lst, key=lambda a: abs(int(a[1]) - mapped_start) + abs(int(a[2]) - mapped_end), reverse=True)
            except Exception:
                print(strand_flag)
                print(ref_lst)
                sys.exit()

            # chr1    585988  827796  ENSG00000230021 .       -       RP5-857K21.4    ensembl_havana  lincRNA
            out_tmp_lst = [":".join([a[-1], a[3], a[6]]) for a in tmp_lst]
            tag_anno_dict[strand_flag] = out_tmp_lst
        for tag in tag_ids:
            new_anno_tag_detail[tag] = tag_anno_dict
        region_anno_detail[":".join(mapped_info)] = tag_anno_dict
    return (map_to_genome_tags, new_anno_tag_detail, region_anno_detail, rebuilt_unanno_info, un_annotated_tags)


def mapped_genome_anno(configure, config_env, mapped_nc_tag_dict, un_annotated_tags_nc, ref_exp, tag_count_dict):
    sample_name = configure.sample
    bowtie_path = config_env["bowtie"]
    tag_fa = configure.tag_fa
    tag_fa_dict = {}
    for i in open(tag_fa, "r").readlines():
        if i.startswith(">"):
            line = i.strip().split("\t")
            tag_fa_dict[line[0].strip(">")] = line[1]
    out_genome_sam = "{0}.genome.sam".format(sample_name)
    out_genome_bam = "{0}.genome.bam".format(sample_name)
    out_genome_sort = "{0}.genome.sort".format(sample_name)
    out_genome_bowtie_stat = "{0}.genome.bowtie.stat".format(sample_name)
    cmd = [
        bowtie_path,
        "-x",
        config_env["genome_index"],
        "-p",
        config_env["cpu_number"],
        config_env["bowtie_para_4_genome"],
        "-f",
        tag_fa,
        "-S",
        out_genome_sam,
        "2>",
        out_genome_bowtie_stat,
        "&&",
        config_env["samtools"],
        "view -bS",
        out_genome_sam,
        ">",
        out_genome_bam,
        "&&",
        config_env["samtools"],
        "sort",
        out_genome_bam,
        "-o",
        out_genome_sort,
        "&&",
        config_env["bedtools"],
        "bamtobed -i",
        out_genome_sort + ".bam",
        ">",
        out_genome_sort + ".bed",
    ]
    cmd = " ".join(cmd)
    print(cmd)
    print("xiegy-----------------------------------map genome")
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    log_message = "align to {0} completed: {1}".format("genome", datetime.datetime.now().strftime("%c"))
    printlg(log_message)
    out_genome_sort_bed_count = out_genome_sort + ".bed.count"
    out_genome_sort_bed_count_handle = open(out_genome_sort_bed_count, "w")
    with open(out_genome_sort + ".bed", "r") as bf:
        for j in bf:
            line = j.strip().split("\t")
            tag_name = line[3]
            tag_count = tag_fa_dict[tag_name]
            line[4] = str(tag_count)
            out_genome_sort_bed_count_handle.write("\t".join(line) + "\n")
    out_genome_sort_bed_count_handle.close()
    cmd = [
        config_env["bedtools"],
        "merge -s -d 10 -c 6,4,5,5 -o first,collapse,sum,median -i",
        out_genome_sort_bed_count,
        ">",
        "{0}.genome.merge.bed".format(sample_name),
    ]
    cmd = " ".join(cmd)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    mapped_anno_bed = "{0}.genome.annotation.info".format(sample_name)
    mapped_unanno_bed = "{0}.genome.unanno.info".format(sample_name)
    mapped_exon_anno_bed = "{0}.exon.annotation.info".format(sample_name)
    cmd = [
        config_env["bedtools"],
        "intersect",
        "-a",
        "{0}.genome.merge.bed".format(sample_name),
        "-b",
        config_env["trans_bed_annotation"],
        "-wa",
        "-wb",
        ">",
        mapped_anno_bed,
        "&&",
        config_env["bedtools"],
        "intersect",
        "-a",
        "{0}.genome.merge.bed".format(sample_name),
        "-b",
        config_env["trans_bed_annotation"],
        "-wa",
        "-v",
        ">",
        mapped_unanno_bed,
        "&&",
        config_env["bedtools"],
        "intersect",
        "-a",
        "{0}.genome.merge.bed".format(sample_name),
        "-b",
        config_env["exon_bed_annotation"],
        "-wa",
        "-wb",
        ">",
        mapped_exon_anno_bed,
    ]
    cmd = " ".join(cmd)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (map_to_genome_tags, new_anno_tag_detail, region_anno_detail, rebuilt_unanno_info, un_annotated_tags) = deal_mapped_info(mapped_anno_bed, mapped_unanno_bed, mapped_nc_tag_dict, un_annotated_tags_nc, ref_exp, tag_count_dict)
    #    with open()
    return (map_to_genome_tags, new_anno_tag_detail, region_anno_detail, rebuilt_unanno_info, un_annotated_tags)


def get_ncRNAs_exp(configure, config_env, ncRNA_lst):
    mapped_nc_tag_dict = {}
    tag_ref_detail = {}
    ref_tag_detail = {}
    ref_exp = {}
    miR_mapped_detail = {}
    mapped_ncRNA_counts = {}
    sample_name = configure.sample
    tag_fa = configure.tag_fa
    distance_re = re.compile(r"NM:i:(\d)")
    mismatch_pattern = re.compile(r"MD:Z:(\d+)[A-Z]?(\d*)")
    bowtie_path = config_env["bowtie"]
    (mature_miRNA, hairpin_info) = deal_mir_info(config_env["mirbase"])
    for n, i in enumerate(ncRNA_lst):
        tag_ref_detail[i] = {}
        ref_tag_detail[i] = {}
        locals()["{0}_sam".format(i)] = "{0}.{1}.sam".format(sample_name, i)
        locals()["{0}_bowtie_stat".format(i)] = "{0}.{1}.bowtie.stat".format(sample_name, i)
        cmd = [
            bowtie_path,
            "-x",
            config_env[i + "_index"],
            "-p",
            config_env["cpu_number"],
            config_env["bowtie_para_4_{0}".format(i)],
            "-f",
            tag_fa,
            "-S",
            locals()["{0}_sam".format(i)],
            "2>",
            locals()["{0}_bowtie_stat".format(i)],
        ]
        cmd = " ".join(cmd)
        print("xiegy----------------------------------get-nc-exp")
        print(cmd)
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        log_message = "align to {0} database completed: {1}".format(i, datetime.datetime.now().strftime("%c"))
        printlg(log_message)
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
    tag_count_dict = get_fa4sRNA(tag_fa)
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
            mir_exp_dict = get_true_miRexp(tag_ref_detail[i], ref_tag_detail[i], mature_miRNA, tag_count_dict)
            ref_exp[i] = mir_exp_dict
    return (ref_exp, mapped_ncRNA_counts, tag_count_dict, mapped_nc_tag_dict, ref_tag_detail)


def main(configure):
    global log_file
    log_file = open(configure.log_file, "w")
    if get_tag_flag(configure.stat_file):
        pass
    config_env = get_config_env(configure.config)
    ncRNA_lst = ["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "scRNA"]
    function_ncRNA_lst = ["miRNA", "piRNA", "snoRNA", "snRNA", "scRNA"]
    (ref_exp, mapped_ncRNA_counts, tag_count_dict, mapped_nc_tag_dict, ref_tag_detail) = get_ncRNAs_exp(configure, config_env, ncRNA_lst)
    total_counts = sum(tag_count_dict.values())
    base_name = os.path.join(configure.outdir, configure.sample)
    un_annotated_tags_nc = {tag_id: tag_count_dict[tag_id] for tag_id in tag_count_dict if tag_id not in mapped_nc_tag_dict}
    (map_to_genome_tags, new_anno_tag_detail, region_anno_detail, rebuilt_unanno_info, un_annotated_tags) = mapped_genome_anno(configure, config_env, mapped_nc_tag_dict, un_annotated_tags_nc, ref_exp, tag_count_dict)
    un_mapped_tag = set(tag_count_dict.keys()) - map_to_genome_tags
    ncRNA_exp_stat = base_name + ".stat"
    tag_mapped_catagory = base_name + ".tag.ncRNA.classification"
    tag_mapped_catagory_handle = open(tag_mapped_catagory, "w")
    tag_mapped_catagory_handle.write("TagId\tTagCount\tCategory\tmapped_item\tregion\n")
    nes = open(ncRNA_exp_stat, "w")
    total_mapped_tags_sum = sum([tag_count_dict[j] for j in map_to_genome_tags])
    total_mapped_nc_tags_sum = sum([tag_count_dict[j] for j in mapped_nc_tag_dict])
    function_nc_tags_sum_ori = sum([mapped_ncRNA_counts.get(j, 0) for j in function_ncRNA_lst])
    total_mapped_tags_ratio = total_mapped_tags_sum / total_counts * 100
    nes.write("#Total tags: {0}\n".format(total_counts))
    nes.write("#Total mapped tags: {0}({1:.2f}%)\n".format(total_mapped_tags_sum, total_mapped_tags_ratio))
    nes.write("#Unmapped tags: {0}({1:.2f}%)\n".format(total_counts - total_mapped_tags_sum, 100 - total_mapped_tags_ratio))
    nes.write("#Category\tMappedTag\tRatio\n")
    exp_cal_category = config_env["RPM"]
    ncRNA_type_list = {"total": ncRNA_lst, "ncRNAall": ncRNA_lst, "func": function_ncRNA_lst, "ncRNAtype": ncRNA_lst}
    exp_cal_method = {"total": total_mapped_tags_sum, "ncRNAall": total_mapped_nc_tags_sum, "func": function_nc_tags_sum_ori}
    for n, i in enumerate(ncRNA_type_list[exp_cal_category]):
        ncRNA_exp_file = base_name + "." + i + ".exp"
        nef = open(ncRNA_exp_file, "w")
        nef.write("GeneSymbol\tTagCount\tRPM\n")
        nc_counts_sum = sum([ref_exp[i][j] for j in ref_exp[i] if j != "other"])
        if exp_cal_category == "ncRNAtype":
            exp_cal_method["ncRNAtype"] = nc_counts_sum

        # ref_tag_detail[i].setdefault(mapped_ref,set()).add(tag_name)
        for j in ref_exp[i]:
            tmp_nc_counts = 0
            if j == "other":
                continue
            if j in ref_tag_detail[i]:
                tmp_mapped_tags = ref_tag_detail[i][j] if j in ref_tag_detail[i] else []
                for t in tmp_mapped_tags:
                    tmp_tag_count = tag_count_dict[t]
                    tag_mapped_catagory_handle.write(t + "\t" + str(tmp_tag_count) + "\t" + i + "\t" + j + "\treference\n")
            tmp_nc_counts = ref_exp[i][j]
            nc_exp = tmp_nc_counts / exp_cal_method[exp_cal_category] * 1000000
            nef.write("{0}\t{1:d}\t{2:.2f}\n".format(j, tmp_nc_counts, nc_exp))
        nef.close()
        nes.write("{0}\t{1:d}\t{2:.2f}%\n".format(i, nc_counts_sum, nc_counts_sum / total_mapped_tags_sum * 100))
        log_message = "calcultion for expression of {0} completed: {1}".format(i, datetime.datetime.now().strftime("%c"))
        printlg(log_message)
    tag_mapped_genome_catagory = base_name + ".tag.genome.classification"
    tag_mapped_genome_catagory_handle = open(tag_mapped_genome_catagory, "w")
    tag_mapped_genome_catagory_handle.write("TagId\tTagCount\tgenebody\topposite\n")
    for tag in new_anno_tag_detail:
        tag_count = tag_count_dict[tag]
        tag_anno_info = new_anno_tag_detail[tag]
        tag_mapped_genome_catagory_handle.write(tag + "\t" + str(tag_count))
        for strand in ["genebody", "opposite"]:
            if strand in tag_anno_info:
                out_mapped_items = ";".join(tag_anno_info[strand])
            else:
                out_mapped_items = "-"
            tag_mapped_genome_catagory_handle.write("\t" + out_mapped_items)
        tag_mapped_genome_catagory_handle.write("\n")
    log_message = "ncRNA expression analysis completed: {0}".format(datetime.datetime.now().strftime("%c"))
    printlg(log_message)
    region_mapped_genome_catagory = base_name + ".region.anno.genome.classification"
    region_mapped_genome_catagory_handle = open(region_mapped_genome_catagory, "w")
    region_mapped_genome_catagory_handle.write("GenomeRegion\tGenebody\topposite\n")
    for region in region_anno_detail:
        region_anno = region_anno_detail[region]
        region_mapped_genome_catagory_handle.write(region)
        for strand in ["genebody", "opposite"]:
            region_mamped_info = ";".join(region_anno[strand]) if strand in region_anno else "-"
            region_mapped_genome_catagory_handle.write("\t" + region_mamped_info)
        region_mapped_genome_catagory_handle.write("\n")
    region_unanno_file = base_name + ".region.unanno.genome.classification"
    region_unanno_file_handle = open(region_unanno_file, "w")
    region_unanno_file_handle.write("Chromosome\tstart\tend\tstrand\tTags\tTotalCount\tMean\n")
    for chr_region in rebuilt_unanno_info:
        region_unanno_file_handle.write(chr_region)
    tag_unmapped = base_name + ".tag.unmapped"
    tag_unmapped_handle = open(tag_unmapped, "w")
    with open(configure.tag_fa, "r") as tmp_f:
        out_flag = 0
        for tmp_line in tmp_f:
            if tmp_line.startswith(">"):
                tag_line = tmp_line.strip(">\n").split("\t")
                if tag_line[0] in un_mapped_tag:
                    tag_unmapped_handle.write(tmp_line)
                    out_flag = 1
            elif out_flag:
                tag_unmapped_handle.write(tmp_line)
                out_flag = 0
    log_message = "Genome annotation completed: {0}".format(datetime.datetime.now().strftime("%c"))
    printlg(log_message)
    nes.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="get annotation for tags, including rRNA, tRNA, miRNA, piRNA, snoRNA, snRNA, scRNA and genome annotation",
    )
    parser.add_argument("-i", "--sample", help="sample name")
    parser.add_argument("-t", "--tag_fa", help="tag file path")
    parser.add_argument("-o", "--outdir", help="outdir path")
    parser.add_argument("-c", "--config", help="configure file path")
    parser.add_argument("-l", "--log_file", help="log file path")
    parser.add_argument("-s", "--stat_file", help="fq stat file path")
    configure = parser.parse_args()
    main(configure)
