#!/usr/bin/env python
# -*- conding:utf-8 -*-
# @AUTHOR: Gui-Yan Xie
# @CONTACT: xieguiyan at hust dot edu dot cn
# @DATE: 2021-11-19 11:28:17
# @DESCRIPTION:
# FILE: plot.py

import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from pathlib import Path

sns.set()


class Plot:
    def __init__(self, inputfile: Path, outputdir: Path):
        self.inputfile = Path(inputfile)
        self.outputdir = outputdir
        self.samprefix = f"{self.outputdir}/{self.inputfile.stem}"

    def read_length_distribution(self) -> None:
        # Reads length distribution
        plt.figure(dpi=300, figsize=(8, 4))
        read = pd.read_table(f"{self.samprefix}.freq.stat", sep="\t", header=None, skipfooter=1)
        read.columns = ["Read length", "count", "Read count percentage"]
        read_len = sns.lineplot(x="Read length", y="Read count percentage", data=read)
        # save image
        read_len.get_figure().savefig(f"{self.outputdir}/distribution_of_read_len.png")
        read_len.get_figure().savefig(f"{self.outputdir}/distribution_of_read_len.pdf")

    def read_type_distribution(self) -> None:
        # Read type distribution
        plt.figure(dpi=300, figsize=(8, 4))
        ncrna_type = pd.read_table(f"{self.samprefix}.stat", header=None, sep="\t", skiprows=4, dtype=str)
        ncrna_type.columns = ["Category", "MappingTag", "Ratio"]
        ncrna_type["Ratio"] = ncrna_type["Ratio"].apply(lambda x: np.nan if x in ["-"] else x[:-1]).astype(float) / 100
        read_type = sns.barplot(x="Category", y="Ratio", data=ncrna_type)
        read_type.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=1))
        # save image
        read_type.get_figure().savefig(f"{self.outputdir}/distribution_of_ncRNA_type.png")
        read_type.get_figure().savefig(f"{self.outputdir}/distribution_of_ncRNA_type.pdf")

    def generate_plot(self):
        self.read_length_distribution()
        self.read_type_distribution()
