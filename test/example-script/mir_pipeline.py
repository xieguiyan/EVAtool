#!/usr/bin/env python
"""
@@ script: sra_proceduere_bowtie.py
@@ description: based on python3.5.2
@@ now based on python3.6.5
@@ author: zhangq
@@ modified by Gui-Yan Xie
"""
import argparse
import datetime
import os
import subprocess
import sys
import tempfile
import re


def println(command, script):
    print(command, file=script)


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


def get_fastq_quals(filename, config_env):
    quals_flag = ""
    solid_re = re.compile(r"T[0-9.]{10,}")
    dump_cmd = [config_env["fastqdump"], "-I", "-X", "1000", "-Z", "--split-spot", filename]
    temp_dump = tempfile.NamedTemporaryFile(delete=False)
    subprocess.run(dump_cmd, stdout=temp_dump, stderr=subprocess.PIPE)
    print(" ".join(dump_cmd))
    temp_dump.close()
    temp_dump_lst = open(temp_dump.name, "r").readlines()
    fastq_lst = temp_dump_lst[1::4]
    solid_flag_num = sum([1 for j in fastq_lst if solid_re.search(j)])
    if solid_flag_num > 100:
        quals_flag = "--dumpbase"
    os.unlink(temp_dump.name)
    print(quals_flag)
    return quals_flag


def generate_command(script, srafiles, workdir, sample, config_file, config_env):
    os.path.isdir(workdir) or os.makedirs(workdir)
    println("cd {workdir}".format(workdir=workdir), script)
    tempdir = tempfile.mkdtemp(
        prefix="analysis.bowtie.", suffix=".{}".format(datetime.datetime.now().strftime("%Y-%m-%d")), dir=workdir
    )
    filter_params = "ILLUMINACLIP:{adpter}:{filter_params}".format(
        adpter=config_env["adp_path"], filter_params=config_env["trimmomatic_sRNA_para"]
    )
    println("mkdir -p {temp} && cd {temp}".format(temp=tempdir), script)
    println("echo 'start processing' > benchmarks.log && date >> benchmarks.log", script)
    println("echo 'start fastq-dump' >> benchmarks.log && date >> benchmarks.log", script)
    sra_lst = []
    for srafile in srafiles.strip().split(","):
        sra = os.path.splitext(os.path.basename(srafile))[0]
        sra_reads = "{sra}_1.fastq.gz".format(sra=sra)
        quals = get_fastq_quals(srafile, config_env)
        println(
            "{fastqdump} {quals} --gzip --split-files {srafile} 1>{sra}.dump.log 2>&1".format(
                fastqdump=config_env["fastqdump"], quals=quals, srafile=srafile, sra=sra
            ),
            script,
        )
        sra_filter_reads_1 = "{sra}.fastq.filter.1.gz".format(sra=sra)
        sra_filter_reads_2 = "{sra}.fastq.filter.fianl.gz".format(sra=sra)
        sra_lst.append(sra_filter_reads_2)
        println(
            "java -jar -Xms8000m -Xmx8000m {trimmomatic} SE -threads {cpu_number} {sra_reads} {out_reads} {filter_params} -trimlog {out_reads}.log 1>{out_reads}_run.log 2>&1".format(
                trimmomatic=config_env["trimmomatic"],
                cpu_number=config_env["cpu_number"],
                sra_reads=sra_reads,
                out_reads=sra_filter_reads_1,
                filter_params=filter_params,
            ),
            script,
        )
        println(
            "java -jar -Xms8000m -Xmx8000m {trimmomatic} SE -threads {cpu_number} {sra_reads} {out_reads} {filter_params} -trimlog {out_reads}.log 1>{out_reads}_run.log 2>&1".format(
                trimmomatic=config_env["trimmomatic"],
                cpu_number=config_env["cpu_number"],
                sra_reads=sra_filter_reads_1,
                out_reads=sra_filter_reads_2,
                filter_params=filter_params,
            ),
            script,
        )
    else:
        println("echo 'end fastq-dump' >> benchmarks.log && date >> benchmarks.log", script)
        println("echo 'start reads2tag translation' >> benchmarks.log && date >> benchmarks.log", script)
        println(
            "{python2} {fq2tag} -i {sras} -c {tag_cut} -o {sample}.fa && echo 'reads2tag translation completed ' >> benchmarks.log && date >> benchmarks.log".format(
                python2=config_env["python2"],
                fq2tag=config_env["fq2tag"],
                sras=",".join(sra_lst),
                tag_cut=config_env["tag_cut"],
                sample=sample,
            ),
            script,
        )
        tag_fa_file = "{0}.fa".format(sample)
        println("echo 'start tag annotation and expression calcultion ' >> benchmarks.log && date >> benchmarks.log", script)
        println(
            "{python3} {tag_anno_exp} -i {sample} -t {tag_fa} -o {tempdir} -c {config_file} -l {sample}.anno.exp.log -s {tag_fa}.freq.stat && echo 'tag annotation completed ' >> benchmarks.log && date >> benchmarks.log".format(
                python3=config_env["python3"],
                tag_anno_exp=config_env["tag_anno_exp"],
                tempdir=tempdir,
                sample=sample,
                tag_fa=tag_fa_file,
                config_file=config_file,
            ),
            script,
        )
    #    'java -jar -Xms8000m -Xmx8000m /project/zhangq/tools/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 SRR5030262_1.fastq SRR5030262_1.filter.11.gz ILLUMINACLIP:/project/zhangq/tools/Trimmomatic-0.36/adapters/sRNA.fa:2:10:4:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 -trimlog SRR5030262_1.filter.log'
    println("echo 'finish processing' >> benchmarks.log && date >> benchmarks.log", script)
    #    println("mv *.txt *.sorted.bam *.log {workdir}".format(workdir=workdir, sra=sra))
    println("cd {workdir}".format(workdir=workdir), script)


#    println("#rm -rf {temp}".format(temp=tempdir))


def main(args):
    script = open(os.path.abspath(args.script), "w")
    config_env = get_config_env(os.path.abspath(args.config))
    generate_command(
        script=script,
        srafiles=args.srafile,
        workdir=os.path.abspath(args.workdir),
        sample=args.sample,
        config_file=args.config,
        config_env=config_env,
    )
    script.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="application description")
    parser.add_argument("-i", "--srafile", help="sra file path, multiple sra should be seperated by comma")
    parser.add_argument("-n", "--sample", help="sample name")
    parser.add_argument("-s", "--script", help="script file path")
    parser.add_argument("-c", "--config", help="configure file path")
    parser.add_argument("-d", "--workdir", help="work directory")
    args = parser.parse_args()
    main(args)
