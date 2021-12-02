#!/usr/bin/env python
#-*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-11-18 21:41:27
# @DESCRIPTION:


#-i: sra file with path (required)
#-o: output directory (required)
#-c: configure file with path ( not required)
#-n: ncRNA type list (not required)

# /home/xiegy/github/EVAtool/test/example-script/example_evatool.sh

/home/xiegy/github/EVAtool/venv/bin/python \
  /home/xiegy/github/EVAtool/evatool/main.py \
  -i /home/xiegy/github/EVAtool/test/example-data/SRR10078125.sra \
  -o /home/xiegy/github/EVAtool/test/tmp_result/SRR10078125
  # -c /home/xiegy/github/EVAtool/evatool/resource/configure.json \
  # -n "miRNA" "rRNA" "tRNA" "piRNA" "snoRNA" "snRNA" "YRNA"