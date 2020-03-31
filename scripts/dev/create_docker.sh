#!/bin/bash
#
# GRIDSS docker container
cd ../..
version=$(grep GRIDSS_VERSION= Dockerfile | cut -b 20-)
echo Building GRIDSS $version
docker build --tag gridss/gridss:latest .
docker build --tag gridss/gridss:$version .
docker push gridss/gridss:latest
docker push gridss/gridss:$version
