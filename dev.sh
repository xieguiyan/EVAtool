#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2021-12-10 15:30:12
# @DESCRIPTION:

# Number of input parameters
# python3 -m pip install --upgrade build
# python3 -m pip install --upgrade twine


# python3 -m build
source venv/bin/activate

#renew requirements
pip freeze > ./requirements-prod.txt

rm -rf dist build
python setup.py sdist bdist_wheel
python3 -m twine upload --repository testpypi dist/*

# python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps evatool
# python3 -m evatool
