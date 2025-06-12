#!/bin/bash

source veclim-data-server/.env

docker build --build-arg VEC_PORTE=${VEC_PORTE} \
             -t veclim-data-server-test \
             .
docker run -d \
           -p ${VEC_PORTE}:${VEC_PORT} \
           --restart always \
           --name veclim-data-server-test \
           -v ${VEC_DATA}:${DIR_DATA} \
           veclim-data-server-test