import sys

sys.path.append("../evatool")


def test_fastq():
    """Test fastq.py"""
    from evatool.utils.fastq import Fastq
    from evatool.utils.config import Config
    from evatool.utils.logger import Logger

    fastq = Fastq(inputfile="../../../test/example-data/SRR8185773.sra", outputdir="../../../test/tmp_result", config=Config(), log=Logger())
    # fastq.process_fastq()
