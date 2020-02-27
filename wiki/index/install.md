## hammurabi X installation instruction

Installing hammurabi X can be done manually in two technical steps, or automatically by using the docker image.
The manual installation requires a computing environment preparation step, 
where [GSL](https://www.gnu.org/software/gsl/), 
[FFTW3](http://www.fftw.org/) and [Google Test](https://github.com/google/googletest) 
are required for a full installation, 
for users who do not want to be bothered with debugging or further modification/development, 
the minimal external libraries are GSL and FFTW3.

> We are currently using FFTW3 version 3.3.8. 
But we are catching with the latest update of all dependencies as much as possible, 
so for a more specific version instruction please read also the 
[Travis CI configuration](https://github.com/hammurabi-dev/hammurabiX/blob/master/.travis.yml).

### package dependencies

The GSL is usually considered as one of the most important standard scientific computing libraries.
In hammurabi X we need GSL library mainly for computing special functions, e.g., the synchrotron and gamma functions.
The FFTW library is used for carrying out fast Fourier transforms in generating random fields according to their power spectra.
The Google Test library is used for testing and debugging the package which is not required for pure users, 
while for programmers it needs no further introduction.

### Installation with Source

First of all, make sure the external packages have been installed properly, 
if Google Test has been ignored then no debugging support will be activated.
The installation setup is written in CMake style, 
so please find the CMake lists file in the root directory of hammurabi X and the important lines are:

```
#-------------- customized zone --------------#

SET(CMAKE_CXX_COMPILER "g++")
OPTION(ENABLE_TESTING "Enable testing for this project" ON)
OPTION(ENABLE_TIMING "Enable timing for this porect" ON)
OPTION(ON_DOCKER "Build on docker image" ON)
OPTION(BUILD_SHARED_LIBS "Build shared library" ON)
OPTION(ENABLE_REPORT "Enable verbose report" ON)

#-------------- instruction ------------------#

# ENABLE_TESTING by default ON, for package build testing,
# ENABLE_TIMING is for checking performance,
# ON_DOCKER by defauly ON, 
# we highly recommend Docker image for non-HPC tasks,
# but if install manually,
# switch it off and modify the path hints LISTED BELOW,
# BUILD_SHARED_LIB by default ON,
# you will be overwhelmed by ENABLE_REPORT,
# switch it off for non-testing tasks,
# 
# you have to specify your local paths of external libraries just below here,
# in some special cases you have to modify FIND_PATH/FIND_LIBRARY functions,
#
# in some special cases, you may want to manually fix LFLAGS or CFLAGS,
#
# if you add new modules/derived classes beyond original code,
# please manually add source file paths to SET(SRC_FILES ...) function,
#
# the last resort would be manually calling homemade Makefile building system,
# you can find cached building files in "cache",
#
# we use Google Test for tests,
# Google Test package is assembled INTO testing modules manually,
# you can either install GoogleTest and cp src dir into install path,
# or just download GoogleTest and specify root dir to GTEST_HINTS
#---------------------------------------------#
```

Noitce that the installation path is registered as INSTALL_ROOT_DIR in CMakeLists.txt.
A simple CMaking process can be:
```
$ cd <package root>
$ mkdir build
$ cd build
$ cmake .. && make -j<n> && make install
```

Once installed, you can run the executable directly with:
```
$ hamx [XML parameter file path]
```
to be specific, the hammurabi executable requires an XML parameter file as input.

### python wrapper

The python wrapper can be found at "<root path>/hampyx".
You can use it as a module without installation like how we use it in the [tutorials](https://github.com/hammurabi-dev/hammurabiX/tree/master/tutorials).
Or, you can install it by
```
$ python setup.py install
```

If you encounter "MKL FATAL ERROR: Cannot load neither libmkl_avx.so nor libmkl_def.so", try:
```
export LD_PRELOAD="/path/to/libmkl_core.so:/path/to/libmkl_sequential.so"
```
This problem may happen in python-C++ routine.

### installation with docker image

The docker image is convenient for making fast and light-weighted productions with hammurabiX,
you can build directly with Dockerfile provided along with the source.
