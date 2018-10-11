# Publishes the GRIDSS docker
version=2.0.1
cd ../docker
docker build --tag gridss/gridss:$version .
docker build --tag gridss/gridss:latest .
docker push gridss/gridss:latest
docker push gridss/gridss:$version