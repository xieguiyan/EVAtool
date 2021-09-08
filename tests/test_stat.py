import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_stat():
    from evatool.utils.stat import Stat
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger

    stat_result = Stat(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    stat_result.pocess_stat()


# test_stat()
