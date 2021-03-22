#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-03-30 17:19:02
# @DESCRIPTION:

#-i: sra file path
#-n: sample name
#-s: script file path
#-c: configure file path
#-d: work directory


/home/xiegy/tools/anaconda3/bin/python \
  /home/xiegy/github/EVAtool/test/example-script/mir_pipeline_zhangq_20180902.py \
  -i /home/xiegy/github/EVAtool/test/example-data/SRR10078125.sra \
  -n SRR10078125 \
  -s /home/xiegy/github/EVAtool/test/example-data/SRR10078125.sh \
  -c /home/xiegy/github/EVAtool/refs/config.miRNA.s10.txt \
  -d /home/xiegy/github/EVAtool/test/example-data/SRR10078125

#bash /home/xiegy/github/EVAtool/test/example-data/DRR006758.sh
