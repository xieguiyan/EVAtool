from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="evatool",
    version="0.0.1",
    keyworkds="evatool",
    description="Extracellular vesicle small RNA processing package",
    author="Dr. Xie and Dr. Liu",
    author_email="xieguiyanathustdotedudotcn",
    url="https://github.com/xieguiyan/EVAtool",
    packages=find_packages(where="src"),
    install_requires=["numpy"],
    license="MIT",
    platforms="any",
    package_dir={"evatool": "evatool"},
    python_requires=">=3.5",
)
