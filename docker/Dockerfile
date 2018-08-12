################## BASE IMAGE ######################

FROM ubuntu:16.04

################## METADATA ######################
LABEL base.image="biocontainers:latest" 
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="1.8.1"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"


RUN bash -c 'echo -e "deb http://archive.ubuntu.com/ubuntu/ xenial main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-updates main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-backports main restricted universe multiverse\n\
deb http://archive.ubuntu.com/ubuntu/ xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list'
	
RUN apt-get clean all && \
	apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y \
	bwa \
	openjdk-8-jre-headless \
	r-base \
	wget

RUN useradd -ms /bin/bash gridss

USER gridss

WORKDIR /home/gridss/

RUN wget https://github.com/PapenfussLab/gridss/releases/download/v1.8.1/gridss-1.8.1-gridss-jar-with-dependencies.jar

WORKDIR /data/

ENTRYPOINT ["java", "-ea", "-Xmx16g", "-Dsamjdk.create_index=true", "-Dsamjdk.use_async_io_read_samtools=true", "-Dsamjdk.use_async_io_write_samtools=true", "-Dsamjdk.use_async_io_write_tribble=true", "-Dgridss.gridss.output_to_temp_file=true", "-cp", "/home/gridss/gridss-1.8.1-gridss-jar-with-dependencies.jar", "gridss.CallVariants", "WORKER_THREADS=4"]