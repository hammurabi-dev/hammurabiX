# README

Hammurabi X (reforged Hammurabi in C++11 std)

Re-developed by:

* Jiaxin Wang (SISSA)
* Tess R. Jaffe (NASA)
* Theo Steininger (MPA)
* Torsten A. Ensslin (MPA)

It simulates 

* synchrotron emission
* Faraday depth

with following physical quantities modeled:

* galactic structures
* galactic magnetic fields
* free(thermal) electron fields
* cosmic ray electron fields

Our code contains (or is benefit from) following packages:

* Hammurabi
* YMW16
* GARFIELDS
* TinyXML2
* PyMultiNest

Major improvements:

* in c++11 std
* XML as parameter file style
* parameters handled collectively
* memory and in/outputs handled collectively
* apply Simpson's rule in line-of-sight integration
* YMW16 regular free electron density field
* divergence-free anisotropic random magnetic field generator
* random free electron density field generator (need further development)
* python wrapper/interface to DRAGON code (under development)
* fully parallelized with MPI (under development)

Dependencies:

* Healpix
* FFTW3
* CFITSIO
* GSL

Compiling and running:
```
make -f [makefile_name]
cd bin/
./hamx [paramfile_name]
```
