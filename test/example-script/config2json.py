#!/usr/bin/env python
"""
@@ script: sra_proceduere_bowtie.py
@@ description: based on python3.6.5
@@ author: Gui-Yan Xie
"""
import argparse
import os
import json


def get_config_env(config_file, workdir):
    config_env = {}
    with open(config_file) as cf:
        for i in cf:
            if i.startswith("#"):
                continue
            if i.find("=") > 1:
                fields = [j.strip() for j in i.strip().split("=")]
                config_env[fields[0]] = fields[1]
    config_json = json.dumps(config_env, indent=4)
    config_file_name = "{workdir}/config_transfered.json".format(workdir=workdir)
    with open(config_file_name, "w") as json_file:
        json_file.write(config_json)
    return config_env


def main(args):
    get_config_env(os.path.abspath(args.config), args.workdir)
    print("Done! Your file has been transfered!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description="application description")
    parser.add_argument("-c", "--config", help="configure file path")
    parser.add_argument("-d", "--workdir", help="work directory")
    args = parser.parse_args()
    main(args)
