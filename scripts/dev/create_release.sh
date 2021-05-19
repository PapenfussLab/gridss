#!/bin/bash
set -o errexit -o pipefail -o noclobber -o nounset
cd ../..
version=$(git tag --points-at HEAD | cut -b 2-)
if [[ "${1-missing}" == "-f" ]] ; then
	version=$(git describe --tags --abbrev=0 --match "v[0-9]*" | cut -b 2-)
fi
if [[ "$version" == "" ]] ; then
	echo "current commit is not a tagged GRIDSS release" 2>&1
	exit 1
fi
echo Building GRIDSS $version 2>&1
rm -rf release/
mkdir release
docker build --build-arg GRIDSS_VERSION=$version --target gridss_export_build_artefacts --output type=local,dest=release --progress=plain . && \
docker build --build-arg GRIDSS_VERSION=$version --target gridss -t gridss:$version -t gridss:latest . && \
docker build --build-arg GRIDSS_VERSION=$version --target gridss_minimal -t gridss_minimal:$version -t gridss_minimal:latest .
cd release
tar -czvf gridss-$version.tar.gz *
# TODO: update conda?
echo docker push gridss/gridss:$version gridss/gridss:latest gridss/gridss_minimal:$version gridss/gridss_minimal:latest 



