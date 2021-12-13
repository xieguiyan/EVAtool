FROM rocker/rstudio:latest
LABEL maintainer="Chun-Jie Liu <chunjie.sam.liu@gmail.com>"

RUN apt-get update && apt-get install -y \
  apt-transport-https \
  ca-certificates \
  curl \
  gnupg-agent \
  software-properties-common \
  cmake \
  curl \
  cython \
  g++ \
  gfortran \
  git \
  libblas-dev \
  libgsl-dev \
  liblapack-dev \
  python-dev \
  python-numpy \
  r-base \
  r-cran-nloptr \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip && pip install --upgrade evatool

ENTRYPOINT [ "evatool" ]