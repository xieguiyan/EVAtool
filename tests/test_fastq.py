import sys
from pathlib import Path

sys.path.append("../EVAtool")


def test_fastq():
    """Test fastq.py"""
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger

    fastq = Fastq(
        inputfile="/home/xiegy/github/EVAtool/test/example-data/SRR10078125.sra",
        outputdir="/home/xiegy/github/EVAtool/test/tmp_result/SRR10078125",
        config=Config(configfile="/home/xiegy/github/EVAtool/refs/reference_config.json"),
        log=Logger(),
        ncrna_lst=["miRNA", "rRNA", "tRNA", "piRNA", "snoRNA", "snRNA", "YRNA"],
    )
    fastq.process_fastq()


test_fastq()
