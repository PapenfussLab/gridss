#!/bin/sh

# Setup script

# Download the GRIDSS jar file and dependencies from:
wget https://github.com/PapenfussLab/gridss/releases/download/v1.8.1/gridss-1.8.1-gridss-jar-with-dependencies.jar

# Download and decompress the ENCODE blacklist file:
mkdir data
wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
gunzip ENCFF001TDO.bed.gz
mv ENCFF001TDO.bed data/ENCODE_blacklist_hg19

mkdir output
