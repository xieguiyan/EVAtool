import sys
import os

sys.path.append("../evatool")


def test_config():
    """Test the config file."""

    from evatool.utils.config import Config

    config = Config()
    print(config.config)


test_config()
