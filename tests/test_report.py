import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_report():
    """Test report.y"""
    from evatool.utils.report import Report
    from evatool.utils.logger import Logger
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config

    report_html = Report(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    report_html.prepare_html()


test_report()
