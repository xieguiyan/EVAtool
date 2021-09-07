from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    version="0.0.1",
    name="evatool",
    keyworkds="evatool",
    description=long_description,
    author="Dr. Xie and Dr. Liu",
    author_email="xieguiyanathustdotedudotcn",
    url="https://github.com/xieguiyan/EVAtool",
    packages=["evatools", "evatools.resource"],
    include_package_data=True,
    package_data={"evatools.resource": ["*.json"]},
    install_requires=["numpy", "pandas"],
    license="MIT",
    platforms="any",
    package_dir={"evatool": "evatool"},
    python_requires=">=3.5",
    scripts=["bin/evatool"],
)
