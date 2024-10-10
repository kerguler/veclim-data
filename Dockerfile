# start from base
FROM ubuntu:22.04

LABEL maintainer="Kamil Erguler <kerguler@gmail.com>"

ARG VEC_PORTE

# install system-wide deps for python and node
RUN apt-get -yqq update
RUN apt-get -yqq install python3-pip python3-dev curl gnupg
RUN apt-get -yqq install proj-bin
RUN apt-get -yqq install libgdal-dev

# copy our application code
COPY ./ /opt/VEClim
WORKDIR /opt/VEClim

# fetch app specific deps
RUN pip3 install --no-cache-dir Cython
RUN pip3 install --no-cache-dir -r ./requirements.txt

# run python package fixes
RUN bash perform_fixes.sh

# expose port
EXPOSE ${VEC_PORTE}

# start app
CMD [ "python3", "./veclim-data-server/veclim-python-server.py" ]
