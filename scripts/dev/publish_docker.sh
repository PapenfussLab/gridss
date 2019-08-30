# Publishes the GRIDSS docker
version=2.4.0
cd ../docker
docker build --tag gridss/gridss:$version .
docker build --tag gridss/gridss:latest .
docker push gridss/gridss:latest
docker push gridss/gridss:$version
