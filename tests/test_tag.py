import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_tag():
    from evatool.utils.tag import Tag
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger

    stat_result = Tag(fastq=Fastq(inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR8185773.sra", outputdir="/home/xiegy/github/EVAtool/test/tmp_result", config=Config(), log=Logger()))
    stat_result.pocess_stat()


test_tag()
