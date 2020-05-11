FROM ubuntu:18.04
LABEL base.image="ubuntu:18.04"

RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic main restricted" > /etc/apt/sources.list
RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic-security main restricted" >> /etc/apt/sources.list
RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic-updates main restricted" >> /etc/apt/sources.list
RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic universe multiverse" >> /etc/apt/sources.list
RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic-security universe multiverse" >> /etc/apt/sources.list
RUN echo "deb http://mirror.aarnet.edu.au/ubuntu/ bionic-updates universe multiverse" >> /etc/apt/sources.list

# CRAN ubuntu package repository for the latest version of R
RUN apt-get update ; apt-get install -y apt-transport-https software-properties-common ; apt-get clean ; rm -rf /var/lib/apt/lists/*
RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add -
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get update ; DEBIAN_FRONTEND=noninteractive apt-get install -y openjdk-8-jre-headless samtools bwa time build-essential make r-base libssl-dev libcurl4-openssl-dev libxml2-dev ; rm -rf /var/lib/apt/lists/*

# R packages used by GRIDSS and PURPLE
ENV R_INSTALL_STAGED=false
RUN Rscript -e 'options(Ncpus=16L, repos="https://cloud.r-project.org/");install.packages(c("tidyverse", "assertthat", "testthat", "NMF", "randomForest", "stringdist", "stringr", "argparser", "R.cache", "BiocManager", "Rcpp", "blob", "RSQLite"))'
RUN Rscript -e 'options(Ncpus=16L, repos="https://cloud.r-project.org/");BiocManager::install(ask=FALSE, pkgs=c("copynumber", "StructuralVariantAnnotation", "VariantAnnotation", "rtracklayer", "BSgenome", "Rsamtools", "biomaRt", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene"))'

ENV GRIDSS_VERSION=2.8.3
ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="$GRIDSS_VERSION"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"

RUN mkdir /opt/gridss/
COPY target/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /opt/gridss/
COPY scripts/*.sh /opt/gridss/
COPY scripts/*.R /opt/gridss/

WORKDIR /data/
ENTRYPOINT ["/opt/gridss/gridss.sh"]
