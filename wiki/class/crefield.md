## hammurabi X crefield class

The ``Crefield`` class is designed for implementing various galactic cosmic-ray electron field models, including internal analytic descriptions and external numerical data.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/crefield.h)

- [CRE base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/cre/cre.cc)
- [analytic CRE source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/cre/cre_ana.cc)
- [uniform CRE source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/field/cre/cre_unif.cc)

## function list:

- **``read_field``**
```
# input
const hamvec<3,double> & (target spatial position)
const double & (target energy)
const Param * (parameter set)
const Grid_cre * (cosmic-ray electron grid)
# return
double (cosmic-ray electron flux)
```
> It reads out the cosmic-ray electron flux at given position and energy in the galactic centric Cartesian frame, regardless of whether the distribution is defined by internal analytic model or from external numerical data. When the base class object is instantiated, it will return a null field distribution.

- **``write_field``**
```
# input
const hamvec<3,double> & (target spatial position)
const double & (target energy)
const Param * (parameter set)
# return
double (cosmic-ray electron flux)
```
> It reads out the cosmic-ray electron flux at given position and energy in the galactic centric Cartesian frame, from the known analytic distribution.

- **``read_grid``**
```
# input
const hamvec<3,double> & (target spatial position)
const double & (target energy)
const Param * (parameter set)
const Grid_cre * (cosmic-ray electron grid)
# return
double (cosmic-ray electron flux)
```
> It reads out the cosmic-ray electron flux at given position and energy in the galactic centric Cartesian frame. Linear interpolation is used since the flux distribution is recorded in the ``Grid_cre`` mesh.

- **``read_grid_num``**
```
# input
const hamvec<3,double> & (target spatial position)
const std::size_t & (target energy index)
const Param * (parameter set)
const Grid_cre * (cosmic-ray electron grid)
# return
double (cosmic-ray electron flux)
```
> It reads out the cosmic-ray electron flux at given position and energy index in the galactic centric Cartesian frame. Note that the energy index means the index defined by ``Grid_cre``,
and it reads out the data stored in the numerical mesh without spectral interpolation, but still with spatial interpolation.

- **``write_grid``**
```
# input
const Param * (parameter set)
Grid_cre * (cosmic-ray electron grid)
# return
-
```
> It writes the internally defined analytic cosmic-ray electron flux distribution into ``Grid_cre``. Note that the grid/mesh defined by ``Grid_cre`` is by default linear Cartesian in the spatial domain, but logarithmic in the spectral domain.

- **``flux_norm``**
```
# input
const hamvec<3,double> & (target spatial position)
const Param * (parameter set)
# return
double (cosmic-ray electron flux normalization)
```
> It calculates the cosmic-ray electron flux normalization at given spatial position. This is mainly used in arpproximating the cosmic-ray synchrotron emissivity from analytic modellings.

- **``flux_idx``**
```
# input
const hamvec<3,double> & (target spatial position)
const Param * (parameter set)
# return
double (cosmic-ray electron flux spectral index)
```
> It calculates the cosmic-ray electron flux spectral index at given spatial position (in the galactic centric Cartesian frame), by assuming the spectral index is independent of the energy. This is mainly used in approximating the cosmic-ray synchrotron emissivity from analytic modellings. It is not available for models with energy dependent spectral index, in which case, we can use ``write_grid`` first and then approximate the synchrotron emissivity with ``read_grid`` or ``read_grid_num``.

- **``spatial_profile``**
```
# input
const hamvec<3,double> & (target spatial position)
const Param * (parameter set)
# return
double (cosmic-ray electron flux spatial profile)
```
> It calculates the cosmic-ray electron flux spatial profile at given spatial position (in the galactic centric Cartesian frame). This is mainly used in approximating the cosmic-ray synchrotron emissivity from analytic modellings.

## method list:

For details of calculating the cosmic-ray electron synchrotron emission with either analytic models or numerical descriptions, please turn to the appendix A of:

[hammurabi X: Simulating Galactic Synchrotron Emission with Random Magnetic Fields](https://arxiv.org/abs/1907.00207)
