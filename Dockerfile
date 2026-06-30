FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y --no-install-recommends apt-utils && apt install -y --no-install-recommends \
  ca-certificates \
  git \
  build-essential \
  clang-14 \
  llvm-14 \
  python3 \
  python3-setuptools \
  python3-pip \
  libgsl-dev \
  libmpfr-dev \
  libmpfrc++-dev \
  libomp-14-dev \
  libboost-all-dev \
  libalglib-dev \
  libopenblas-dev \
  && ln -sf /usr/lib/x86_64-linux-gnu/libgsl.so /usr/local/lib/libgsl.so \
  && pip3 install wheel mpmath numpy scipy seaborn pandas matplotlib \
  && pip3 install --index-url https://download.pytorch.org/whl/cpu torch

RUN mkdir -p /OFP-Repair
COPY . /OFP-Repair
WORKDIR /OFP-Repair
