#!/bin/bash


mount_dir=$(pwd)

chmod a+x "${mount_dir}"/docker/build_me.sh

cd docker && ./build_me.sh && cd ../

parentdir="$(dirname "$mount_dir")"

docker run \
  --name scrna_pipeline \
  --user $(id -u):$(id -g) \
  -d \
  --rm \
  -v "${parentdir}":/pipeline \
  ubuntur35:pipeline
