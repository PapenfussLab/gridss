################## BASE IMAGE ######################

FROM ubuntu:18.04

################## METADATA ######################
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="2.4.0"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"


RUN bash -c 'echo -e "deb http://archive.ubuntu.com/ubuntu/ xenial main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-updates main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-backports main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list'

RUN apt-get clean all && apt-get update
RUN apt-get install -y openjdk-8-jre-headless
RUN apt-get install -y r-base
RUN apt-get install -y bwa
RUN apt-get install -y time

RUN useradd -ms /bin/bash gridss
RUN mkdir /data
RUN chown gridss /data
RUN chgrp users /data
RUN chmod 777 /data

RUN mkdir /data/gridss
ADD https://github.com/PapenfussLab/gridss/releases/download/v2.4.0/gridss-2.4.0-gridss-jar-with-dependencies.jar /data/gridss/
ADD https://github.com/PapenfussLab/gridss/releases/download/v2.4.0/gridss.sh /data/gridss/
ADD https://github.com/PapenfussLab/gridss/releases/download/v2.4.0/gridss_lite.sh /data/gridss/
RUN chmod +x  /data/gridss/gridss*.sh
WORKDIR /data/gridss/

ENTRYPOINT ["/data/gridss/gridss.sh", "--jar", "/data/gridss/gridss-2.4.0-gridss-jar-with-dependencies.jar"]
