import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_stat():
    """Test stat.py"""
    from evatool.utils.stat import Stat
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger
    from evatool.utils.fastq import Fastq
    from evatool.utils.tag import Tag

    stat = Stat(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    stat.stat_match()


test_stat()
