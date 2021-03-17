FROM ubuntu:20.04
LABEL base.image="ubuntu:20.04"

# Setup CRAN ubuntu package repository
RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y \
		wget \
		apt-transport-https \
		software-properties-common \
		dirmngr \
		gnupg && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*
RUN apt-key adv \
	--keyserver hkp://keyserver.ubuntu.com:80 \
	--recv-keys 0xE298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y \
		apt-utils \
		gawk \
		openjdk-11-jre-headless \
		samtools \
		bwa \
		bedtools \
		r-base \
		time \
		build-essential \
		autotools-dev \
		autoconf \
		autogen \
		make \
		libssl-dev \
		libcurl4-openssl-dev \
		libxml2-dev \
		zlib1g-dev \
		libbz2-dev \
		liblzma-dev && \
	apt-get clean && \
	rm -rf /var/lib/apt/lists/*

# R packages used by GRIDSS
ENV R_INSTALL_STAGED=false
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/");install.packages(c( "tidyverse", "assertthat", "testthat", "randomForest", "stringdist", "stringr", "argparser", "R.cache", "BiocManager", "Rcpp", "blob", "RSQLite" ))'
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/");BiocManager::install(ask=FALSE, pkgs=c( "copynumber", "StructuralVariantAnnotation", "VariantAnnotation", "rtracklayer", "BSgenome", "Rsamtools", "biomaRt", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene" ))'

ENV GRIDSS_VERSION=2.11.0
ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="$GRIDSS_VERSION"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"

RUN mkdir /opt/gridss/
#COPY target/github_package/* /opt/gridss/
# downloading from github allows the docker images to be built with zero reliance on the local environment
RUN wget -O /opt/gridss/gridss-${GRIDSS_VERSION}.tar.gz https://github.com/PapenfussLab/gridss/releases/download/v${GRIDSS_VERSION}/gridss-${GRIDSS_VERSION}.tar.gz
RUN tar zxvf /opt/gridss/gridss-${GRIDSS_VERSION}.tar.gz -C /opt/gridss/ && \
	tar zxvf /opt/gridss/gridsstools.src.tar.gz -C /opt/gridss/ && \
	cd /opt/gridss/src/main/c/gridsstools/htslib && \
	autoconf && autoheader && make -j 8 && \
	cd /opt/gridss/src/main/c/gridsstools/ && \
	autoconf && autoheader && make  -j 8 gridsstools && \
	cp gridsstools /opt/gridss/gridsstools
RUN chmod +x /opt/gridss/gridsstools /opt/gridss/*.R /opt/gridss/*.sh
ENV PATH="/opt/gridss:$PATH"

WORKDIR /data/

