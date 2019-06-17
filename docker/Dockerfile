################## BASE IMAGE ######################

FROM ubuntu:16.04

################## METADATA ######################
LABEL base.image="biocontainers:latest"
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
RUN apt-get install -y wget

RUN useradd -ms /bin/bash gridss
RUN mkdir /data
RUN chown gridss /data
RUN chgrp users /data
RUN chmod 777 /data

USER gridss
RUN mkdir /data/gridss
WORKDIR /data/gridss
RUN wget https://github.com/PapenfussLab/gridss/releases/download/v2.4.0/gridss-2.4.0-gridss-jar-with-dependencies.jar

# TODO: check java version in container.
# Must be at least java 8 191
# -XX:ActiveProcessorCount
# -XX:MaxRAMPercentage
ENTRYPOINT ["java",
	"-XX:+UnlockExperimentalVMOptions", 
	"-XX:+UseCGroupMemoryLimitForHeap", 
	"-XX:MaxRAMFraction=1",
	"-XshowSettings:vm",
	"-Dsamjdk.create_index=true",
	"-Dsamjdk.use_async_io_read_samtools=true",
	"-Dsamjdk.use_async_io_write_samtools=true",
	"-Dsamjdk.use_async_io_write_tribble=true",
	"-Dgridss.gridss.output_to_temp_file=true",
	"-Dsamjdk.buffer_size=4194304",
	"-cp",
	"/data/gridss/gridss-2.4.0-gridss-jar-with-dependencies.jar",
	"gridss.CallVariants"]
