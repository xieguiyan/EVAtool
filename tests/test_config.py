import sys
import os


sys.path.append("../EVAtool")


def test_config():
    from evatool.utils.config import Config

    config = Config()
    print(config.config)


test_config()
