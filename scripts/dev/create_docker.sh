#!/bin/bash
#
# GRIDSS docker container
cd ../..
pom_version=$(grep "gridss</version>" pom.xml | cut -f 2 -d '>' | cut -f 1 -d '-')
version=$(grep GRIDSS_VERSION= Dockerfile | cut -b 20-)
if [[ $pom_version != $version ]] ; then
	echo Dockerfile GRIDSS version does not match pom.xml
	exit 1
fi
echo Building GRIDSS $version
docker build --tag gridss/gridss:latest .
docker build --tag gridss/gridss:$version .
docker push gridss/gridss:latest
docker push gridss/gridss:$version
