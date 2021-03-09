#!/bin/bash
conda deactivate
set -o errexit -o pipefail -o noclobber -o nounset
cd ../../
version=$(grep "gridss</version>" pom.xml | cut -f 2 -d '>' | cut -f 1 -d '-')
echo Packaging GRIDSS $version
rm -rf target
mvn clean package
mkdir target/github_package
cp scripts/*.sh target/github_package/
cp scripts/*.R target/github_package/
cp target/gridss-$version-gridss-jar-with-dependencies.jar target/github_package/
cd src/main/c/gridsstools/htslib
autoheader
autoconf
./configure && make clean && make -j $(nproc)
cd -
cd src/main/c/gridsstools
autoheader
autoconf
./configure && make clean && make -j $(nproc) all
cd -
cp src/main/c/gridsstools/gridsstools target/github_package
cd src/main/c/gridsstools/htslib
make clean
cd -
cd src/main/c/gridsstools
make clean
cd - 
tar -czvf target/github_package/gridsstools.src.tar.gz src/main/c/*
cd target/github_package
tar -czvf gridss-$version.tar.gz *
