#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-09-03 15:29:17
# @DESCRIPTION:


class Config(object):
    def __init__(self, config_file):
        self.config_file = config_file
        self.config = {}
        self.read_config()

    def read_config(self):
        with open(self.config_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if line == "":
                    continue
                key, value = line.split("=")
                self.config[key] = value

    def get(self, key):
        return self.config[key]
