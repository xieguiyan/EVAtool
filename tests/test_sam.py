import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_sam():
    """Test sam.py"""
    from evatool.utils.sam import SAM
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger
    from evatool.utils.fastq import Fastq

    sam = SAM(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    sam.ncRNA_map()
