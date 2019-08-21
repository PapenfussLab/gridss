################## BASE IMAGE ######################

FROM ubuntu:18.04

################## METADATA ######################

ENV GRIDSS_VERSION 2.5.1
ENV SAMBAMBA_VERSION 0.7.0

ENV GRIDSS_JAR "/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar"

LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="${GRIDSS_VERSION}"
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

RUN useradd -Ms /bin/bash gridss

# Install sambamba
RUN mkdir -p /opt/sambamba
ADD https://github.com/biod/sambamba/releases/download/v${SAMBAMBA_VERSION}/sambamba-${SAMBAMBA_VERSION}-linux-static.gz /opt/sambamba/

RUN cd /opt/sambamba \
    && gunzip sambamba-${SAMBAMBA_VERSION}-linux-static.gz \
    && mv sambamba-${SAMBAMBA_VERSION}-linux-static sambamba \
    && chmod +x sambamba

# Install gridss to /opt/gridss, and hard link it from /data/ to maintain backwards compatability
RUN mkdir -p /opt/gridss
ADD https://github.com/PapenfussLab/gridss/releases/download/v${GRIDSS_VERSION}/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar ${GRIDSS_JAR}
ADD https://github.com/PapenfussLab/gridss/releases/download/v${GRIDSS_VERSION}/gridss.sh /opt/gridss/
# ADD https://github.com/PapenfussLab/gridss/releases/download/v${GRIDSS_VERSION}/gridss_lite.sh /opt/gridss/

RUN chmod +x  /opt/gridss/gridss*.sh

RUN chown -R gridss /opt
RUN chgrp users /opt
RUN chmod -Rv 755 /opt/gridss/

RUN mkdir -p /data/gridss
RUN ln /opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /data/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
RUN ln /opt/gridss/gridss.sh /data/gridss/gridss.sh
# RUN ln /opt/gridss/gridss-lite.sh /data/gridss/gridss-lite.sh

# RUN chmod -R 777 /data

WORKDIR /opt/gridss/
ENV PATH="/opt/gridss:/opt/sambamba:${PATH}"

# See https://stackoverflow.com/a/37904830/2860731 for "echo $GRIDSS_JAR" explanation
ENTRYPOINT ["/opt/gridss/gridss.sh"]
