#!/bin/bash
set -o errexit -o pipefail -o noclobber -o nounset
cd ../../
version=$(grep "gridss</version>" pom.xml | cut -f 2 -d '>' | cut -f 1 -d '-')
echo Packaging GRIDSS $version
rm -rf target/github_package
mvn clean package
mkdir target/github_package
cp scripts/*.sh target/github_package/
cp scripts/*.R target/github_package/
cp target/gridss-$version-gridss-jar-with-dependencies.jar target/github_package/
cd -
cd src/main/c/gridsstools
./configure && make clean && make
cd -
cp src/main/c/gridsstools/gridsstools target/github_package
cd target/github_package
tar -czvf gridss-$version.tar.gz *
