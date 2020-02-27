## hammurabi X bfield class

The ``Bfield`` class is designed for implementing various galactic magnetic field models,
including both the regular(mean) and random(turbulent) components.
In fact, the base class ``Bfield`` itself is just a cover, where the two derived ``Breg`` and ``Brnd`` are the true base classes which defines some meaningful functions.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/bfield.h)

- [regular magnetic field base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/breg.cc)
- [uniform regular magnetic field source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/breg_unif.cc)
- [LSA model (regular magnetic field) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/breg_lsa.cc)
- [Jaffe model (regular magnetic field) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/breg_jaffe.cc)
- [random magnetic field base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/brnd.cc)
- [ES model (random magnetic field) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/brnd_es.cc)
- [parametric MHD model (random magnetic field) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/b/brnd_mhd.cc)

## function list:

### for regular fields

- **``Breg::read_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_breg * (regular magnetic field grid)
# return
hamvec<3,double> (magnetic field vector)
```
> It reads out the regular magnetic field vector at the given galactic centric Cartesian frame position, the returned field vector is described in the same Cartesian frame. It is flexible to handle either built-in analytic field model (call ``write_field``) or external numerical input data (call ``read_grid``).

- **``Breg::write_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
# return
hamvec<3,double> (magnetic field vector)
```
> It writes out the regular magnetic field vector at given galactic centric Cartesian frame position with built-in analytic field model. This function offers magnetic field vector directly without bypassing any internal storage nor interpolation.

- **``Breg::read_grid``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_breg * (regular magnetic field grid)
# return
hamvec<3,double> (magnetic field vector)
```
> It reads out the regular magnetic field vector at given galactic centric Cartesian frame position with external numerical input data, which has been imported from disk and stored inside the ``Grid_breg`` object. Linear interpolation of a 3D vector field is used since the ``Grid_breg`` defines only a finite number of sample points within a Cartesian mesh.

- **``Breg::write_grid``**
```
# input
const Param * (parameter set)
Grid_breg * (regular magnetic field grid)
# return
-
```
> It records the built-in analytic regular magnetic field distribution to the ``Grid_breg`` sotrage, which can be either exported to an external file or used internally.

### for random fields

- **``Brnd::read_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_brnd * (random magnetic field grid)
# return
hamvec<3,double> (magnetic field vector)
```
> It reads out the random magnetic field vector at the given galactic centric Cartesian frame position, the returned field vector is described in the same Cartesian frame. It only handles either built-in or external numerical input data (call ``read_grid``). When the base-class object is instantiated, it returns a null vector.

- **``Brnd::read_grid``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_brnd * (random magnetic field grid)
# return
hamvec<3,double> (magnetic field vector)
```
> It reads out the random magnetic field vector at given galactic centric Cartesian frame position, from the ``Grid_brnd`` regardless the origin of the data stored inside.

- **``Brnd::write_grid``**
```
# input
const Param *
const Breg *
const Grid_breg *
Grid_brnd *
# return
-
```
> It is the function that implement the internal random field generators. Once the random magnetic field is generated according to the Cartesian grid defined by ``Grid_brnd``, it will be stored and later can be used by ``read_grid``. The ``Breg`` and ``Grid_breg`` are listed as input argument in case the random component generation depends on the regular magnetic field distribution.

## method list:

For the general algorithms of random magnetic field generators, please turn to:

[hammurabi X: Simulating Galactic Synchrotron Emission with Random Magnetic Fields](https://arxiv.org/abs/1907.00207)

Some captions can be found along with the source file, and we will provde more comprehensive descriptions upon requests.
