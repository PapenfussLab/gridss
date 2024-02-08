ARG UBUNTU_VERSION=20.04
FROM ubuntu:$UBUNTU_VERSION AS gridss_base_closest_mirror
# Use the closest mirror so apt-get doesnt take ages
RUN sed -i -e 's/http:\/\/archive\.ubuntu\.com\/ubuntu\//mirror:\/\/mirrors\.ubuntu\.com\/mirrors\.txt/' /etc/apt/sources.list

# Set up a C build environment for gridsstools, samtools, and R packages
FROM gridss_base_closest_mirror AS gridss_c_build_environment
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y ca-certificates
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
	libssl-dev \
	libcurl4-openssl-dev \
	libxml2-dev \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libdeflate-dev \
	build-essential \
	autotools-dev \
	autoconf \
	autogen \
	make \
	wget \
	libomp-dev \
    python3-pip \
	&& rm -rf /var/lib/apt/lists/* \
	&& pip3 install pysam

# compile gridsstools
FROM gridss_c_build_environment AS gridss_builder_c
RUN mkdir /opt/gridss/
ARG GRIDSS_VERSION
COPY src/main/c /opt/gridss/src/main/c
COPY src/test/resources /opt/gridss/src/test/resources
RUN cd /opt/gridss/src/main/c/gridsstools/htslib && \
	autoreconf -i && ./configure && make -j 8 && \
	cd .. && \
	autoreconf -i && ./configure && make -j 8 && \
	cp gridsstools /opt/gridss/

# compile GRIDSS Java code
FROM maven:3.8.4-jdk-11 AS gridss_builder_java
RUN mkdir /opt/gridss/
WORKDIR /opt/gridss/
# Download maven dependencies first so docker can cache them
COPY pom.xml /opt/gridss/
COPY repo /opt/gridss/repo
# run all stages so all dependencies are cached
RUN mvn -Dmaven.artifact.threads=8 verify && rm -rf target
# Build GRIDSS jar
ARG GRIDSS_VERSION
COPY src /opt/gridss/src
RUN mvn -T 1C -Drevision=${GRIDSS_VERSION} package && \
	cp target/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /opt/gridss/

FROM gridss_c_build_environment AS gridss
# Setup CRAN ubuntu package repository
# apt-get clean not required for ubuntu images
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
		apt-transport-https \
		software-properties-common \
		dirmngr \
		gnupg && \
	apt-key adv \
	--keyserver hkp://keyserver.ubuntu.com:80 \
	--recv-keys 0xE298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
	apt-get update && apt-get install -y \
		apt-utils \
		gawk \
		openjdk-11-jre-headless \
		bwa \
		hmmer \
		bedtools \
		bcftools \
		r-base \
		time \
		libomp-dev \
		perl-modules \
		libtext-soundex-perl \
		python3-h5py \
		rsync \
		curl \
		libxml2-dev \
		libcairo2-dev \
		libgit2-dev \
		default-libmysqlclient-dev \
		libpq-dev \
		libsasl2-dev \
		libsqlite3-dev \
		libssh2-1-dev \
		libxtst6 \
		libharfbuzz-dev \
		libfribidi-dev \
		libfreetype6-dev \
		libpng-dev \
		libtiff5-dev \
		libjpeg-dev \
		unixodbc-dev \
	&& rm -rf /var/lib/apt/lists/*
# samtools needs to be installed from source since the OS package verion is too old
RUN mkdir /opt/samtools && \
	cd /opt/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
	tar -jxf samtools-1.14.tar.bz2 && \
	cd samtools-1.14 && \
	autoheader && \
	autoconf -Wno-syntax && \
	./configure && \
	make install && \
	cd ~ && \
	rm -rf /opt/samtools
### Repeat Masker and dependencies
RUN mkdir /opt/trf && \
	cd /opt/trf && \
	wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 && \
	chmod +x trf*.linux64 && \
	ln -s trf*.linux64 trf
# Turns out we need makeblastdb as well as rmblastn (https://github.com/PapenfussLab/gridss/issues/535)
RUN mkdir /opt/rmblast && \
	cd /opt/rmblast && \
	wget https://www.repeatmasker.org/rmblast/rmblast-2.14.1+-x64-linux.tar.gz && \
	tar --no-anchored --strip-components 2 -xvzf rmblast-2.14.1+-x64-linux.tar.gz rmblastn makeblastdb && \
	rm rmblast-2.14.1+-x64-linux.tar.gz
RUN cd /opt/ && \
	wget http://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz && \
	tar zxf RepeatMasker-*.tar.gz && \
	rm RepeatMasker-*.tar.gz
# Install GATK
RUN mkdir /opt/gatk && \
	cd /opt/gatk && \
	 wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && \
      unzip gatk-4.2.6.1.zip && \
      rm gatk-4.2.6.1.zip
RUN apt update && apt --yes install default-jdk
### Kraken2 and dependencies
# dustmasker from e-direct: (or is this in ncbi-blast as well?)
RUN mkdir /opt/blast && \
	cd /opt/blast && \
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz && \
	tar zxf ncbi-blast-*.tar.gz && \
	mv ncbi-blast-*/bin/* . && \
	rm -r ncbi-blast-*
ENV KRAKEN_VERSION=2.1.2
RUN mkdir /opt/kraken2 && \
	cd /opt/kraken2 && \
	wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v$KRAKEN_VERSION.tar.gz && \
	tar zxf v*.tar.gz && \
	cd kraken2* && \
	./install_kraken2.sh /opt/kraken2 && \
	cd .. && \
	rm -r kraken2-$KRAKEN_VERSION v*.tar.gz
RUN sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" && \
	mv $HOME/edirect /opt/edirect
ENV PATH="/opt/gridss/:/opt/RepeatMasker:/opt/rmblast/:/opt/trf:/opt/kraken2:/opt/blast:/opt/edirect:/opt/gatk:$PATH"
# configure repeatmasker
RUN cd /opt/RepeatMasker && \
	perl configure \
		-default_search_engine rmblast \
		-rmblast_dir /opt/rmblast \
		-trf_prgm /opt/trf/trf \
		-hmmer_dir /usr/local/bin
# R packages used by GRIDSS - R package need the C toolchain installed
ENV R_INSTALL_STAGED=false
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/");install.packages(c( "tidyverse", "assertthat", "testthat", "randomForest", "stringdist", "stringr", "argparser", "R.cache", "BiocManager", "Rcpp", "blob", "RSQLite" ))'
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/");BiocManager::install(ask=FALSE, pkgs=c( "copynumber", "StructuralVariantAnnotation", "VariantAnnotation", "rtracklayer", "BSgenome", "Rsamtools", "biomaRt", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg19" ))'
# Install GRIDSS
ARG GRIDSS_VERSION
ENV GRIDSS_VERSION=${GRIDSS_VERSION}
ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
LABEL software="GRIDSS"
LABEL software.version="$GRIDSS_VERSION"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"
RUN mkdir /opt/gridss/ /data
COPY --from=gridss_builder_c /opt/gridss/gridsstools /opt/gridss/
COPY --from=gridss_builder_java /opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /opt/gridss/
COPY scripts/gridss \
	scripts/gridss_annotate_vcf_kraken2 \
	scripts/gridss_annotate_vcf_repeatmasker \
	scripts/gridss_extract_overlapping_fragments \
	scripts/gridss_somatic_filter \
	scripts/virusbreakend \
	scripts/virusbreakend-build \
	scripts/gridss.config.R \
	scripts/libgridss.R \
	scripts/link_breakpoints \
	scripts/revert_sup_low_mapq_ua_alignment.py \
	/opt/gridss/
RUN chmod +x /opt/gridss/* && \
	chmod -x /opt/gridss/*.R
WORKDIR /data/

# Copy build artifact locally
FROM scratch AS gridss_export_build_artefacts
ARG GRIDSS_VERSION
COPY --from=gridss /opt/gridss/* /

