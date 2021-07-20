#!/bin/bash
conda deactivate
set -o errexit -o pipefail -o noclobber -o nounset
cd ../..
version=$(git tag --points-at HEAD | cut -b 2-)
if [[ "${1-missing}" == "-f" ]] ; then
	version=$(git describe --tags --abbrev=0 --match "v[0-9]*" | cut -b 2-)
fi
if [[ "$version" == "" ]] ; then
	version=$(git symbolic-ref --short HEAD)-$(git rev-parse --short HEAD)
fi

echo Building GRIDSS $version 2>&1
rm -rf release/
mkdir release
(cd src/main/c/gridsstools && make clean)
(cd src/main/c/gridsstools/htslib && make clean)
(cd src/main/c/ && tar czf gridsstools-src-$version.tar.gz gridsstools)
mv src/main/c/gridsstools-src-$version.tar.gz release/
cp LICENSE release/ # https://bioconda.github.io/contributor/linting.html#gpl-requires-license-distributed
docker build --build-arg GRIDSS_VERSION=$version --target gridss_export_build_artefacts --output type=local,dest=release . # --progress=plain
docker build --build-arg GRIDSS_VERSION=$version --target gridss -t gridss/gridss:$version -t gridss/gridss:latest .
docker build --build-arg GRIDSS_VERSION=$version --target gridss_minimal -t gridss/gridss_minimal:$version -t gridss/gridss_minimal:latest .
docker build --build-arg GRIDSS_VERSION=$version --target gridss -t gridss/virusbreakend:$version -t gridss/virusbreakend:latest .
cd release
chmod +x *
chmod -x *.R *.jar
tar czf gridss-$version.tar.gz *
echo docker push gridss/gridss:$version 
echo docker push gridss/gridss:latest
echo docker push gridss/gridss_minimal:$version
echo docker push gridss/gridss_minimal:latest 
echo docker push gridss/virusbreakend:$version
echo docker push gridss/virusbreakend:latest 


