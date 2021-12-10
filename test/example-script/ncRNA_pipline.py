#!/usr/bin/env python
"""
@@ script: sra_proceduere_bowtie.py
@@ description: based on python3.6.5
@@ author: Gui-Yan Xie
"""
import argparse
import datetime
import os
import subprocess
import tempfile
import re
import json
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a file handler
handler = logging.FileHandler("benchmarks.log")
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
# add the file handler to the logger
logger.addHandler(handler)


def generate_command(srafile, workdir, config_file):
    resultdir = "{dir}/analysis.bowtie.{time}".format(dir=workdir, time=datetime.datetime.now().strftime("%Y-%m-%d"))
    isExist = os.path.isdir(resultdir)
    if not isExist:
        os.makedirs(resultdir)
    else:
        pass
    os.chdir(resultdir)
    with open(config_file, "r") as c:
        config_env = json.load(c)
    logger.info("Directory has been created!\n")
    sra = os.path.splitext(os.path.basename(srafile))[0]
    sra_reads = "{sra}_1.fastq.gz".format(sra=sra)
    quals = "--dumpbase"
    transfq = "{fastqdump} {quals} --gzip --split-files {srafile} 1>{sra}.dump.log 2>&1".format(fastqdump=config_env["fastqdump"], quals=quals, srafile=srafile, sra=sra)
    runfq = subprocess.run(transfq, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
    if runfq.returncode == 0:
        logger.info("sucess in transfer sra to fastq!")
    else:
        logger.info("error in transfer sra to fastq!")
    sra_filter_reads_1 = "{sra}.fastq.filter.1.gz".format(sra=sra)
    filter_params = "ILLUMINACLIP:{adpter}:{filter_params}".format(adpter=config_env["adp_path"], filter_params=config_env["trimmomatic_sRNA_para"])
    trimfq = "java -jar -Xms8000m -Xmx8000m {trimmomatic} SE -threads {cpu_number} {sra_reads} {out_reads} {filter_params} -trimlog {out_reads}.log 1>{out_reads}_run.log 2>&1".format(
        trimmomatic=config_env["trimmomatic"],
        cpu_number=config_env["cpu_number"],
        sra_reads=sra_reads,
        out_reads=sra_filter_reads_1,
        filter_params=filter_params,
    )
    runtrim = subprocess.run(trimfq, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
    if runtrim.returncode == 0:
        logger.info("sucess in trimmomatic for fastq!")
    else:
        logger.info("error in trimmomatic for fastq!")
    logger.info("start reads2tag translation")
    fq2tag = "{python2} {fq2tag} -i {sratrim} -c {tag_cut} -o {sra}.fa".format(
        python2=config_env["python2"],
        fq2tag=config_env["fq2tag"],
        sratrim=sra_filter_reads_1,
        tag_cut=config_env["tag_cut"],
        sra=sra,
    )
    runfq2tag = subprocess.run(fq2tag, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
    if runfq2tag.returncode == 0:
        logger.info("reads2tag translation completed!")
    else:
        logger.info("error in reads2tag translation!")
    logger.info("start tag annotation and expression calcultion!")
    tag_fa_file = "{0}.fa".format(sra)
    taganno = "{python3} {tag_anno_exp} -i {sra} -t {tag_fa} -o {resultdir} -c {config_file} -l {sra}.anno.exp.log -s {tag_fa}.freq.stat".format(
        python3=config_env["python3"],
        tag_anno_exp=config_env["tag_anno_exp"],
        resultdir=resultdir,
        sra=sra,
        tag_fa=tag_fa_file,
        config_file=config_file,
    )
    logger.info("Tag annotation completed!")
    os.chdir(workdir)


def main(args):
    generate_command(srafile=args.srafile, workdir=os.path.abspath(args.workdir), config_file=args.config)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="application description")
    parser.add_argument("-i", "--srafile", help="sra file path, only for one file")
    parser.add_argument("-c", "--config", help="configure file path")
    parser.add_argument("-d", "--workdir", help="work directory")
    args = parser.parse_args()
    main(args)
