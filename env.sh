#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-12-28 14:39:45
# @DESCRIPTION:

# Number of input parameters
git pull

[ -d venv ] || eval '`which python3` -m venv venv'


source venv/bin/activate
pip install --upgrade pip
pip install -r requirements-dev.txt