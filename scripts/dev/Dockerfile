FROM maven:3.6.3-jdk-8

# Update apt
# Install build essentials
# Then htslib prerequisites
RUN apt-get update -y && \
    apt-get install -y \
        build-essential \
        autotools-dev \
        autoconf \
        autogen \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*