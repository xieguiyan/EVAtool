#!/usr/bin/env bash
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan@hust.edu.cn
# @DATE: 2021年 08月 31日 星期二 21:22:28 CST
# @DESCRIPTION:

#-i: sra file path
#-n: sample name
#-s: script file path
#-c: configure file path
#-d: work directory


/home/xiegy/tools/anaconda3/bin/python \
  /home/xiegy/github/EVAtool/test/example-script/ncRNA_pipline.py \
  -i /home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra \
  -n SRR8185773 \
  -s /home/xiegy/github/EVAtool/test/example-data/SRR8185773.sh \
  -c /home/xiegy/github/EVAtool/refs/config_ncRNA.txt \
  -d /home/xiegy/github/EVAtool/test/example-data/SRR8185773
