from io import DEFAULT_BUFFER_SIZE
from setuptools import setup, find_packages


MAJOR = 0
MINOR = 1
MICRO = 1
ISRELEASED = True
VERSION = f"{MAJOR}.{MINOR}.{MICRO}"

DESCRIPTION = "Extracellular Vesicles small RNAs Abundance and quantification tool"


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="evatool",
    version=VERSION,
    keyworkds="evatool",
    description=DESCRIPTION,
    long_description=long_description,
    author="Dr. Xie and Dr. Liu",
    author_email="xieguiyan@hust.eud.cn",
    url="https://github.com/xieguiyan/EVAtool",
    packages=find_packages(),
    include_package_data=True,
    package_data={"evatool.resource": ["*.json"]},
    install_requires=["numpy", "pandas", "seaborn"],
    license="MIT",
    platforms="any",
    package_dir={"evatool": "evatool"},
    python_requires=">=3.5",
)
