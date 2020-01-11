#!/bin/bash

docker build -f Dockerfile.ub -t ubuntu1804:base .
docker build -f Dockerfile.r -t ubuntu:r35 .
docker build -f Dockerfile.pipe -t ubuntur35:pipeline .
