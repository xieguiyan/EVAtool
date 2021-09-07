#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at.hust.edu.cn
# @DATE: 2021-09-03 15:29:17
# @DESCRIPTION:


import json
from pathlib import Path

current_path = Path(__file__).parent


class Config(object):
    def __init__(self, configfile: Path = current_path / "../resource/configure.json"):
        self.configfile = Path(configfile)
        self.config = self.read_config()

    def read_config(self):
        try:
            with open(self.configfile, "r") as f:
                return json.load(f)
        except IOError:
            print(f"{self.configfile} not exists.")


# config = Config().config
