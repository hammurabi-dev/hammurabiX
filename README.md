# README

Next Generation Hammuarbi

Open to HAMMURABI/DRAGON/IMAGINE developing groups ONLY!

Currently, we simulate synchrotron emission, Faraday depth with following physical quantities modeled:

* galactic structures
* galactic magnetic fields
* free(thermal) electron fields
* cosmic ray electron fields

This code contains (or is benefit from) following packages:

* Hammurabi
* YMW16
* GARFIELDS
* TinyXML2
* PyMultiNest

Major improvements:

* in c++11 std
* XML as parameter file style
* parameters handled by Pond class collectively
* memory, in/outputs handled by Grid class collectively
* YMW16 regular free electron density field
* divergence-free random magnetic field generator (need further development)
* random free electron density field generator (need further development)
* not producing dust emission (need further development)
* python wrapper
* simple mcmc pipeline with PyMultiNest for regular fields
* fully parallelized (under development)

Dependencies:

* Healpix (Healpy)
* FFTW3
* CFITSIO
* GSL
* PyMultiNest