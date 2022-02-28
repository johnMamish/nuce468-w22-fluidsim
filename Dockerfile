FROM nvidia/cuda:11.0-devel-ubuntu20.04

# thanks for killing centos, redhat <3
# now we have to fiddle with the mirrorlist to do any package modifications
# should probably just switch to the ubuntu image, but i don't really want to re-download it
# RUN sed -i -e "s|mirrorlist=|#mirrorlist=|g" /etc/yum.repos.d/CentOS-*
# RUN sed -i -e "s|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g" /etc/yum.repos.d/CentOS-*

# RUN dnf install -y git

ENV DEBIAN_FRONTEND noninteractive

RUN apt update
RUN apt install -y git cmake python3

ENV CUDA_INSTALL_PATH /usr/local/cuda-11.0

COPY . /CE468_labs

RUN make all -C /CE468_labs

CMD [ "make", "-C", "/CE468_labs", "test" ]
