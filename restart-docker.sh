#!/bin/bash

source .env

docker build --build-arg VEC_PORTE=${VEC_PORTE} \
             -t ${VEC_NAME} \
             .
docker run -d \
           -p ${VEC_PORTE}:${VEC_PORT} \
           --restart always \
           --name ${VEC_NAME} \
           -v ${VEC_DATA}:${DIR_DATA} \
           ${VEC_NAME}