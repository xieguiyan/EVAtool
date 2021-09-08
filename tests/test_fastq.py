import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_fastq():
    """Test fastq.py"""
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger

    fastq = Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger())
    fastq.process_fastq()


# test_fastq()
