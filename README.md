# hammurabi README
##### (current version X.alpha)

hammurabi is a HEALPix-based tool for simulating observables, 
such as **polarized synchrotron** and **thermal dust emission**, 
from models of physical inputs such as **magnetic fields**, **dust** and **electron distributions**, etc.

It is a modular C++ framework into which you can add your own models easily, 
and then use it to perform the line-of-sight integration to compute the observables. 
Please cite the original [Waelkens et al. (2009)](https://www.aanda.org/articles/aa/abs/2009/08/aa10564-08/aa10564-08.html) paper if you use hammurabi.

#### about X.alpha release
*Currently we are in alpha testing release,
with only minimal amount of features available.
A full release of version X is scheduled in 2018.*


X.alpha provides simulation of: 

* polarized synchrotron emission
* Faraday depth

with following physical modelings:

* galactic structures
* galactic magnetic fields
* free electron fields
* cosmic ray electron fields

with support from following packages:

* [YMW16](https://bitbucket.org/psrsoft/ymw16)
* [GARFIELDS](https://academic.oup.com/mnras/article-lookup/doi/10.1111/j.1365-2966.2008.13341.x)
* [TinyXML2](https://github.com/leethomason/tinyxml2)

with major improvements:

* in c++11 std
* XML style parameter file 
* parameters, memory, and I/O handled collectively
* apply Simpson's rule in integration
* YMW16 regular free-electron field
* divergence-free anisotropic turbulent magnetic field generator
* turbulent free-electron field generator

with package dependencies:

* [Healpix](https://healpix.jpl.nasa.gov/)
* [FFTW3](http://www.fftw.org/)
* [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
* [GSL](https://www.gnu.org/software/gsl/)

#### compiling:
```
$ cd [package root]
$ make -f install/[makefile]
```

#### running:
```
$ cd bin/
$ ./hamx [paramfile]
```
remarks: we suggest users to compile without -DNDEBUG to verify pipeline, while run the code with -NDEBUG to save computing time.

#### contact:
*bug reports and code contributions are warmly welcomed,
feel free to contact Jiaxin Wang, Tess Jaffe, and Torsten Ensslin*
