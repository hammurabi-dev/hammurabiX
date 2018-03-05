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
RUN wget https://sourceforge.net/projects/healpix/files/latest/download
RUN mkdir healpix
RUN tar zxf download -C ./healpix --strip-components=1
WORKDIR healpix
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
ENV HEALPIX /home/lab/healpix

# CLEAN PACKAGES
RUN rm -f *.tar.gz

# HAMMURABI X
RUN git clone https://bitbucket.org/hammurabicode/hamx
WORKDIR hamx
# move to slave branch (if necessary)
RUN git fetch origin unitest
RUN git checkout unitest
RUN make test
RUN make
RUN make install
WORKDIR ..
