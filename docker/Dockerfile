FROM ubuntu:latest

MAINTAINER jiaxin.wang@sissa.it

RUN apt-get update
RUN apt-get install -y build-essential git gfortran autoconf gsl-bin libgsl-dev wget

RUN mkdir /home/lab
WORKDIR /home/lab

# CFITSIO
RUN wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz && tar xzf cfitsio_latest.tar.gz
WORKDIR cfitsio
RUN ./configure --prefix=/usr/local/ && make && make install 
WORKDIR ..

# FFTW
RUN wget http://www.fftw.org/fftw-3.3.7.tar.gz && tar xzf fftw-3.3.7.tar.gz
WORKDIR fftw-3.3.7
RUN ./configure --enable-threads --enable-openmp --enable-shared --prefix=/usr/local/ && make && make install
WORKDIR ..

# HEALPIX
RUN wget http://downloads.sourceforge.net/project/healpix/Healpix_3.31/Healpix_3.31_2016Aug26.tar.gz && tar xzf Healpix_3.31_2016Aug26.tar.gz 
WORKDIR Healpix_3.31
# write config file on fly
RUN echo '4\n\
	/sur/local/bin\n\
	/usr/local/include\n\
	4\n\
	y\n\
	0\n'\
	> hlpx_config
RUN ./configure < hlpx_config && make 
WORKDIR ..
ENV HEALPIX_TARGET optimized_gcc 
ENV HEALPIX /home/Downloads/Healpix_3.31

# CLEAN PACKAGES
RUN rm -f *.tar.gz

# HAMMURABI X
RUN git clone https://bitbucket.org/hammurabicode/hamx
WORKDIR hamx
# on fix branch
RUN git fetch origin fix
RUN git checkout fix
RUN make
WORKDIR ..
