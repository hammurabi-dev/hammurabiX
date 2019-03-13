# hammurabi README
##### (version X.02.02)

hammurabi is an open-source (GNU General Public License v3) HEALPix-based tool 
for simulating full-sky Galactic foreground observables.
Outputs of hammurabi include:

- **polarized synchrotron emission**
- **thermal dust emission** [in progress]
- **dispersion measure** 
- **Faraday depth**

Essential physical inputs/modelings required during simulation include:  

* **Galactic magnetic fields**
* **dust distribution** 
* **(cosmic-ray/thermal) electron distribution**

hammurabi is a modular C++ framework which is friendly to user defined models.
In version X, we mainly focus on improving its numerical reliablility and scalability.

Please check our [**WIKI PAGE**](https://bitbucket.org/hammurabicode/hamx/wiki/Home) for more detailed technical information.

Original hammurabi source code can be found [**here**](https://sourceforge.net/projects/hammurabicode/).

### hammurabi team publications:

- [Simulating polarized Galactic synchrotron emission at all frequencies. The Hammurabi code](https://www.aanda.org/articles/aa/abs/2009/08/aa10564-08/aa10564-08.html)

### contact
*bug reports and code contributions are warmly welcomed, feel free to contact*

- [Tess Jaffe](https://science.gsfc.nasa.gov/sed/bio/tess.jaffe)
- [Torsten Ensslin](https://wwwmpa.mpa-garching.mpg.de/~ensslin/)
- [Jiaxin Wang](http://www.sissa.it/app/members.php?ID=222)

## version info:

- currently in [TESTING](./tests) release with limited amount of features available
- synchrotron emission with radial integration test presented in wiki
- multi-threading performace presented in wiki

### patch X.02.02
- fix integration precision behavior

### patch X.02.01
- observable output resolution can be defined independently for each map
- redesigned ouput map dictionary entry in hampyx
- mem-leak fixed

### update X.02
- in EARLY release
- customisable radial shell thickness
- thermal dust emission not implemented yet
- free-free absorption removed
- independent vector module: [hvec.h](./include/hvec.h)
- (crude) python wrapper: [hampyx](./hampyx)

### notice with hampyx

- If you encounter "MKL FATAL ERROR: Cannot load neither libmkl_avx.so nor libmkl_def.so", try:

```
export LD_PRELOAD="/path/to/libmkl_core.so:/path/to/libmkl_sequential.so"
```

this problem may happen in python-C++ routine
