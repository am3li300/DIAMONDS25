#syntax=docker/dockerfile:1
FROM ubuntu
LABEL AUTHOR="Mert Erden"
RUN apt-get update
RUN apt-get -y install git
RUN apt-get -y install python3-pip
RUN git clone https://github.com/merterden98/ADAGIO
RUN cd ADAGIO && pip install . 
RUN git clone https://github.com/kap-devkota/GLIDER.git
RUN cd GLIDER/glide && pip install .