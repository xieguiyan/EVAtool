#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-11-18 15:36:19
# @DESCRIPTION:

from jinja2 import Environment, FileSystemLoader
import datetime
from pathlib import Path

template_path = Path.cwd() / "evatool" / "resource"


class Report:
    def __init__(self, inputfile: Path, outputdir: Path, config: Path, ncrna_list: list):
        self.inputfile = Path(inputfile)
        self.outputdir = outputdir
        self.config = config
        self.ncrna_list = ncrna_list
        self.samprefix = f"{self.outputdir}/{self.inputfile.stem}"

    def generate_html(self, body, body2, stoptime, img_path, config, sam_info):
        env = Environment(loader=FileSystemLoader(template_path))
        template = env.get_template("template_report.html")
        with open(f"{self.outputdir}/Report_result.html", "w+") as fout:
            html_content = template.render(stop_time=stoptime, body=body, body2=body2, img_path=img_path, config=config, sam_info=sam_info)
            fout.write(html_content)

    def load_readlen_data(self):
        len_dis = {}
        with open(f"{self.samprefix}.freq.stat", "r") as freg:
            for i in freg:
                if i.startswith("ok"):
                    len_stat = i.strip().split("\t")
                    len_dis["fi"] = len_stat[1]
                    len_dis["thr"] = len_stat[2]
                    len_dis["fo"] = len_stat[3]
                    len_dis["fif"] = len_stat[4]
                    len_dis["out_reads"] = len_stat[5]
                    len_dis["ratio"] = len_stat[6]
                    len_dis["total"] = len_stat[7]
        return len_dis

    def load_ncrnatype_data(self):
        map_info = {}
        all_map_info = []
        with open(f"{self.samprefix}.stat", "r") as type_freg:
            for t in type_freg:
                if t.startswith("#"):
                    anno_map_info = t.strip("#|\n").split(":")
                    map_info[anno_map_info[0]] = anno_map_info[1]
                else:
                    ncrna_stat = t.strip().split("\t")
                    if ncrna_stat[0] != "Category":
                        ncrna_dis = {}
                        ncrna_dis["Category"] = ncrna_stat[0]
                        ncrna_dis["MappedTag"] = ncrna_stat[1]
                        ncrna_dis["Ratio"] = ncrna_stat[2]
                        all_map_info.append(ncrna_dis)
        return all_map_info

    def prepare_html(self):
        body = []
        body2 = []
        img_path = {}
        sam_info = {}
        len_dis = self.load_readlen_data()
        all_map_info = self.load_ncrnatype_data()
        img_path["read_len"] = "distribution_of_read_len.png"
        img_path["ncrna_type"] = "distribution_of_ncRNA_type.png"
        time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        body.append(len_dis)
        body2 = all_map_info
        sam_info["sam_name"] = self.inputfile.stem
        sam_info["sam_path"] = self.outputdir
        sam_info["ncrna_list"] = self.ncrna_list
        config = self.config
        self.generate_html(body, body2, time, img_path, config, sam_info)
