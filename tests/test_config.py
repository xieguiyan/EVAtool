import sys
import os


sys.path.append("../EVAtool")


def test_config():
    """Test the config file."""
    print(sys.path)
    import evatool

    from evatool.utils.config import Config

    config = Config()
    print(config.config)


test_config()
