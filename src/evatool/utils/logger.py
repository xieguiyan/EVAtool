#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 15:33:56
# @DESCRIPTION:


class Logger(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_file_handler = open(self.log_file, "w")

    def __del__(self):
        self.log_file_handler.close()

    def log(self, msg):
        self.log_file_handler.write(msg + "\n")
