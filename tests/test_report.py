import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_report():
    """Test report.y"""
    from evatool.utils.report import Report

    report_html = Report(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result")
    report_html.prepare_html()


test_report()
