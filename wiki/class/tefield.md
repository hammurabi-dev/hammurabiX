## hammurabi X tefield class

The ``TEfield`` class is designed for implementing various galactic thermal electron field models,
including both the regular(mean) and random(turbulent) components.
In fact, the base class ``TEfield`` itself is just a cover, where the two derived ``TEreg`` and ``TErnd`` are the true base classes which defines some meaningful functions.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/tefield.h)

- [regular thermal electron base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/te/tereg.cc)
- [YMW16 model (regular thermal electron) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/te/tereg_ymw16.cc)
- [uniform (regular thermal electron) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/te/tereg_unif.cc)
- [random thermal electron base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/te/ternd.cc)
- [default (random theral electron) source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/te/ternd_dft.cc)

## function list:

### for regular fields

- **``TEreg::read_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_tereg * (regular thermal electron field grid)
# return
double (thermal electron field density)
```
> It reads out the regular thermal electron field density at the given galactic centric Cartesian frame position, the returned field density is described in the same Cartesian frame. It is flexible to handle either built-in analytic field model (call ``write_field``) or external numerical input data (call ``read_grid``).

- **``TEreg::write_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
# return
double (thermal electron field density)
```
> It writes out the regular thermal electron field density at given galactic centric Cartesian frame position with built-in analytic field model. This function offers thermal electron field density directly without bypassing any internal storage nor interpolation.

- **``TEreg::read_grid``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_tereg * (regular thermal electron field grid)
# return
double (thermal electron field density)
```
> It reads out the regular thermal electron field density at given galactic centric Cartesian frame position with external numerical input data, which has been imported from disk and stored inside the ``Grid_tereg`` object. Linear interpolation of a 3D density field is used since the ``Grid_tereg`` defines only a finite number of sample points within a Cartesian mesh.

- **``TEreg::write_grid``**
```
# input
const Param * (parameter set)
Grid_tereg * (regular thermal electron field grid)
# return
-
```
> It records the built-in analytic regular thermal electron field distribution to the ``Grid_tereg`` sotrage, which can be either exported to an external file or used internally.

### for random fields

- **``TErnd::read_field``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_ternd * (random thermal electron field grid)
# return
double (thermal electron field density)
```
> It reads out the random thermal electron field density at the given galactic centric Cartesian frame position, the returned field density is described in the same Cartesian frame. It only handles either built-in or external numerical input data (call ``read_grid``). When the base-class object is instantiated, it returns a null density.

- **``TErnd::read_grid``**
```
# input
const hamvec<3,double> & (target position)
const Param * (parameter set)
const Grid_ternd * (random thermal electron field grid)
# return
double (thermal electron field density)
```
> It reads out the random thermal electron field density at given galactic centric Cartesian frame position, from the ``Grid_ternd`` regardless the origin of the data stored inside.

- **``TErnd::write_grid``**
```
# input
const Param *
const TEreg *
const Grid_tereg *
Grid_ternd *
# return
-
```
> It is the function that implement the internal random field generators. Once the random thermal electron field is generated according to the Cartesian grid defined by ``Grid_ternd``, it will be stored and later can be used by ``read_grid``. The ``TEreg`` and ``Grid_tereg`` are listed as input argument in case the random component generation depends on the regular thermal electron field distribution.

## method list:

Currently we haven't put the random thermal electron field into scientific applications, so in the future the explicit implementation will be changed (but the API remains unchanged).
