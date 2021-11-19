import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_plot():
    """Test plot.py"""
    from evatool.utils.plot import Plot

    plot_result = Plot(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result")
    plot_result.generate_plot()


test_plot()
