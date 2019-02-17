# hammurabi README
##### (version X.02)

hammurabi is an open-source (GNU General Public License v3) HEALPix-based tool 
for simulating full-sky Galactic foreground observables.
Outputs of hammurabi include:

* **polarized synchrotron emission**
* **thermal dust emission** [in progress]
* **dispersion measure** 
* **Faraday depth**

Essential physical inputs/modelings required during simulation include:  

* **Galactic magnetic fields**
* **dust distribution** 
* **(cosmic-ray/thermal) electron distribution**

hammurabi is a modular C++ framework which is friendly to user defined models.
In version X, we mainly focus on improving its numerical reliablility and scalability.

Please check our [**WIKI PAGE**](https://bitbucket.org/hammurabicode/hamx/wiki/Home) for more detailed technical information.

Original hammurabi source code can be found [**here**](https://sourceforge.net/projects/hammurabicode/).

Please cite the original [Waelkens et al. (2009)](https://www.aanda.org/articles/aa/abs/2009/08/aa10564-08/aa10564-08.html) paper if you use hammurabi.

### contact
*bug reports and code contributions are warmly welcomed, feel free to contact*

* [Tess Jaffe](https://science.gsfc.nasa.gov/sed/bio/tess.jaffe)
* [Torsten Ensslin](https://wwwmpa.mpa-garching.mpg.de/~ensslin/)
* [Jiaxin Wang](http://www.sissa.it/app/members.php?ID=222)

## version info:

### update X.02
* in EARLY release
* thermal dust emission not implemented yet
* free-free absorption removed
* independent vector module: [hvec.h](./include/hvec.h)
* (crude) python wrapper: [hampyx](./hampyx)

### version X.01
*Currently we are in [TESTING](./tests) release,
with limited amount of features available.
Not all modelings have been thoroughly cross-checked with original hammurabi!*
