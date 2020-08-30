#!/bin/bash
#
cd ../../
version=$(grep "gridss</version>" pom.xml | cut -f 2 -d '>' | cut -f 1 -d '-')
echo Packaging GRIDSS $version
rm -r target/github_package
mvn clean package
mkdir target/github_package
cp scripts/*.sh target/github_package/
cp scripts/*.R target/github_package/
cp target/gridss-$version-gridss-jar-with-dependencies.jar target/github_package/
cd target/github_package
tar -czvf gridss-$version.tar.gz *
