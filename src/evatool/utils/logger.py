#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan.at.gmail.com
# @DATE: 2021-09-03 15:33:56
# @DESCRIPTION:

import logging


class Logger(object):
    def __init__(self, testlog):
        self.testlog = testlog
        self.logtofile = self.mylogger()

    def mylogger(self):
        logger = logging.getLogger("my_logger")
        logger.setLevel(logging.INFO)
        fh = logging.FileHandler("/home/xiegy/github/EVAtool/src/evatool/utils/log.txt")
        fh.setLevel(logging.INFO)
        print(self)

        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        fh.setFormatter(formatter)

        logger.addHandler(fh)

        logger.info(self)


# Logger.mylogger("success in log functions!")
