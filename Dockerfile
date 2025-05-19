FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y

RUN apt-get -y install build-essential

RUN apt-get -y install clang

RUN apt-get -y install libboost-dev

RUN apt-get -y install libomp-dev 

#RUN apt-get -y install python3.8

#RUN apt-get -y install python3-tk

#RUN apt-get -y install python3-pil

#RUN apt-get -y install python3-pil.imagetk

COPY docker/. /usr/src/myapp

WORKDIR /usr/src/myapp

RUN mkdir bin

RUN make clean

RUN make all

COPY data/. /usr/src/myapp/data

CMD ["bin/main", "/usr/src/myapp/data/data.csv", "/usr/src/myapp/data/params.txt", "/usr/src/myapp/out" ]
