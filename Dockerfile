FROM ubuntu:20.04
LABEL base.image="ubuntu:20.04"

# CRAN ubuntu package repository for the latest version of R
RUN apt-get update && \
    apt-get install -y \
        dirmngr \
        gnupg \
        apt-transport-https \
        ca-certificates \
        software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Standard installations
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    apt-utils \
    openjdk-8-jre-headless \
    samtools \
    bwa \
    bedtools \
    time \
    build-essential \
    make \
    r-base \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# R packages used by GRIDSS and PURPLE
ENV R_INSTALL_STAGED=false
RUN Rscript -e 'options(Ncpus=16L, repos="https://cloud.r-project.org/"); install.packages(c("tidyverse", "assertthat", "testthat", "randomForest", "stringdist", "stringr", "argparser", "R.cache", "BiocManager", "Rcpp", "blob", "RSQLite"))'
RUN Rscript -e 'options(Ncpus=16L, repos="https://cloud.r-project.org/"); BiocManager::install(ask=FALSE, pkgs=c("copynumber", "StructuralVariantAnnotation", "VariantAnnotation", "rtracklayer", "BSgenome", "Rsamtools", "biomaRt", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "Biobase"))'
# NMF must be installed AFTER Biobase
RUN Rscript -e 'options(Ncpus=16L, repos="https://cloud.r-project.org/"); install.packages(c("NMF"))'

ENV GRIDSS_VERSION=2.10.2
ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
LABEL version="1"
LABEL software="GRIDSS"
LABEL software.version="$GRIDSS_VERSION"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"

RUN mkdir /opt/gridss/
COPY target/github_package/* /opt/gridss/
RUN chmod +x /opt/gridss/gridsstools /opt/gridss/*.R /opt/gridss/*.sh
ENV PATH="/opt/gridss:$PATH"

WORKDIR /data/
ENTRYPOINT ["/opt/gridss/gridss.sh"]
