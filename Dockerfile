FROM nvidia/cuda:11.0-devel-ubuntu20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt update
RUN apt install -y git cmake python3 python3-opencv python3-pil python3-matplotlib gdb

ENV CUDA_INSTALL_PATH /usr/local/cuda-11.0

COPY . /CE468_labs

RUN make all -C /CE468_labs

CMD [ "make", "-C", "/CE468_labs", "test" ]
