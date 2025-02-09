FROM ubuntu:22.04 as base
LABEL authors="nicholashathaway"

ARG DEBIAN_FRONTEND="noninteractive"
ARG SHELL="/bin/bash"
ARG LANG="en_US.UTF-8"
ARG LANGUAGE="en_US.UTF-8"
ARG LC_ALL="en_US.UTF-8"
ARG CPU_COUNT=15
ARG TIME_ZONE=Etc/UTC
ARG SEEKDEEP_VERSION=develop

RUN ulimit -n 10000

# install basics
RUN apt-get update
RUN apt-get -yq dist-upgrade
RUN apt-get install -yq --no-install-recommends autotools-dev autoconf libtool automake build-essential curl wget file git locales libssl-dev libcurl4-gnutls-dev ca-certificates xz-utils zlib1g-dev libbz2-dev liblzma5 liblzma-doc liblzma-dev openssh-client gcc-10 g++-10 python3-dev python3-pip ccache git-lfs

# install common bio tools
RUN apt-get install -yq --no-install-recommends lastz ncbi-blast+ ncbi-tools-bin hmmer minimap2 bwa bowtie2 cmake samtools mummer muscle pigz

# install snakemake
RUN pip3 install snakemake

# set environment locale
RUN echo "$LANG UTF-8" >> /etc/locale.gen
RUN echo "LANG=$LANG" > /etc/locale.conf
RUN echo "LC_ALL=$LC_ALL" >> /etc/environment
RUN echo "LANGUAGE=$LANGUAGE" >> /etc/environment
RUN locale-gen $LANG
RUN update-locale LANG=$LANG

RUN DEBIAN_FRONTEND=noninteractive TZ=$TIME_ZONE apt-get -y install tzdata

# add github key so can git clone from github
RUN mkdir ~/.ssh/
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts

# install the desired version of github
WORKDIR /opt
RUN git clone https://github.com/bailey-lab/SeekDeep.git
WORKDIR /opt/SeekDeep
RUN git checkout develop
RUN ./setup.py  --libs SeekDeep:${SEEKDEEP_VERSION} --numCores ${CPU_COUNT} --overWrite --symlinkBin

# clean up
RUN rm -fr external/build

# expose ports
EXPOSE 8000-9999


# update path
ENV PATH "/opt/SeekDeep/bin:$PATH"

# reset working directory
WORKDIR /


