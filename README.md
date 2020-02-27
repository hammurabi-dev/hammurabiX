# hammurabi X README

[![Build Status](https://travis-ci.org/hammurabi-dev/hammurabiX.svg?branch=master)](https://travis-ci.org/hammurabi-dev/hammurabiX)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.01889/status.svg)](https://doi.org/10.21105/joss.01889)

hammurabi is an open-source (GNU General Public License v3) package
for simulating full/partial-sky Galactic emissions.
The outputs of (current version) hammurabi X include:

- **polarized synchrotron emission**
- **dispersion measure** 
- **Faraday depth**

Essential physical inputs/modelings required during simulation include:  

- **Galactic magnetic field**
- **(cosmic-ray & thermal) electron distribution**

hammurabi X is a modular C++ framework which is friendly to user defined models.
Based on earlier versions of hammurabi, we mainly focus on improving its numerical reliablility and scalability.

Please check our [**WIKI**](https://github.com/hammurabi-dev/hammurabiX/tree/master/wiki) for more detailed technical information.

[**hammurabi version 3**](https://sourceforge.net/projects/hammurabicode/) is still accessible.

### hammurabi team publications:

- [hammurabi X: a C++ package for simulating Galactic emissions](https://joss.theoj.org/papers/10.21105/joss.01889#)

- [hammurabi X: Simulating Galactic Synchrotron Emission with Random Magnetic Fields](https://iopscience.iop.org/article/10.3847/1538-4365/ab72a2)

- [Simulating polarized Galactic synchrotron emission at all frequencies. The Hammurabi code](https://www.aanda.org/articles/aa/abs/2009/08/aa10564-08/aa10564-08.html)

### contact
*bug reports and code contributions are warmly welcomed, feel free to contact*

- [Jiaxin Wang](https://gioacchinowang.github.io/)
- [Tess Jaffe](https://science.gsfc.nasa.gov/sed/bio/tess.jaffe)
- [Torsten Ensslin](https://wwwmpa.mpa-garching.mpg.de/~ensslin/)

### acknowledgements

- We copied and modified functions/classes from HEALPix.
- We use TinyXML2 as the XML file parser.
