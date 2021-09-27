import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_bam():
    """Test bam.py"""
    from evatool.utils.bam import Bam
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger
    from evatool.utils.fastq import Fastq

    bam = Bam(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    bam.process_bam()


test_bam()
