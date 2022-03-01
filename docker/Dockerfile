#FROM python:3.10
#RUN apt-get update && pip install numpy pandas
FROM ubuntu:latest
RUN apt update \
  && apt -y install software-properties-common dirmngr apt-transport-https lsb-release ca-certificates \
  && add-apt-repository ppa:deadsnakes/ppa \
  && apt-get update \
  && apt-get install -y python3-pip python3 \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip \
  && pip3 install numpy pandas xlrd
#  && update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 1 \
#  && apt-get -y install --reinstall python3-numpy \
#  && pip3 install numpy --upgrade
