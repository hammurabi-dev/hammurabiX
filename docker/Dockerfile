FROM ubuntu:latest

MAINTAINER jiaxin.wang@sjtu.edu.cn

RUN apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --fix-missing build-essential libblas-dev liblapack-dev gsl-bin libgsl-dev cmake
RUN apt-get install -y --fix-missing git wget vim unzip
RUN apt-get install -y --fix-missing python3 python3-dev python3-pip
RUN pip3 install numpy scipy matplotlib astropy healpy jupyter pytest nbval

RUN mkdir /home/lab
WORKDIR /home/lab

# FFTW3
RUN mkdir /home/lab/fftw
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz && tar xzf fftw-3.3.8.tar.gz -C fftw --strip-components 1
WORKDIR /home/lab/fftw
RUN ./configure --enable-threads --enable-openmp --enable-shared --prefix=/tmp/local/ && make && make install
WORKDIR /home/lab

# GOOGLE TEST
RUN mkdir /home/lab/googletest
RUN wget https://github.com/google/googletest/archive/v1.10.x.tar.gz && tar xzf v.1.10.x.tar.gz -C googletest --strip-components 1
WORKDIR /home/lab/googletest
RUN mkdir /home/lab/googletest/build
WORKDIR /home/lab/googletest/build
RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/tmp/local .. && make && make install
RUN mkdir /tmp/local/src/googletest && cp ../googletest/src/* /tmp/local/src/googletest
WORKDIR /home/lab

# HAMMURABI X
RUN git clone https://github.com/hammurabi-dev/hammurabiX.git hamx
WORKDIR /home/lab/hamx
RUN mkdir /home/lab/hamx/build
WORKDIR /home/lab/hamx/build
RUN cmake .. && make install
WORKDIR /home/lab/hamx
RUN python3 setup.py install
WORKDIR /home/lab

# CLEAN PACKAGES
RUN rm -f *.tar.gz

# SET PATH
ENV PATH /tmp/local/hammurabi/bin:${PATH}
ENV LD_LIBRARY_PATH /tmp/local/hammurabi/lib:/tmp/local/lib:${LD_LIBRARY_PATH}
