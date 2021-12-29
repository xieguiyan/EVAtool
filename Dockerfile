FROM rocker/rstudio:latest
ENV DEBIAN_FRONTEND noninteractive
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
  git \
  libblas-dev \
  libgsl-dev \
  liblapack-dev \
  python-dev \
  python-numpy \
  zlib1g-dev \
  && apt-add-repository universe \
  && apt-get update \
  && add-apt-repository ppa:openjdk-r/ppa \
  && apt-get install -y openjdk-8-jdk \
  && apt-get install -y python3-pip \
  && rm -rf /var/lib/apt/lists/*

# java
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

# RUN pip install --upgrade pip && pip install --upgrade evatool
RUN  python3 -m pip install --upgrade pip &&  python3 -m pip install evatool --upgrade

ENTRYPOINT [ "evatool" ]