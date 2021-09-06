#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at.hust.edu.cn
# @DATE: 2021-09-03 15:33:56
# @DESCRIPTION:

import logging
from pathlib import Path


class Logger(object):
    def __init__(self, message, logfile: Path):
        self.message = message
        self.logfile = logfile
        self.logtofile = self.mylogger()

    def mylogger(self):
        logger = logging.getLogger("my_logger")
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler(self.logfile)
        fh.setLevel(logging.INFO)

        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(formatter)

        logger.addHandler(fh)

        logger.info(self.message)


# Logger("Success in log functions!", "/home/xiegy/github/EVAtool/src/evatool/utils/log.txt").logtofile
