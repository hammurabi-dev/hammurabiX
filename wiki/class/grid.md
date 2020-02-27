## hammurabi X Grid class

The ``Grid`` class is designed for manipulating **numerical** fields and observables, and interfacing the memory and disk data collectively.
Data type of each field is by default defined as array/arrays of doubles.
There are multiple source files for hosting the base and derived classes' implementations.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/grid.h)
- [base class source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid.cc)
- [regular magnetic grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_breg.cc)
- [random magnetic grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_brnd.cc)
- [regular thermal electron grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_tereg.cc)
- [random thermal electron grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_ternd.cc)
- [cosmic ray electron grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_cre.cc)
- [regular dust grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_dreg.cc)
- [random dust grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_drnd.cc)
- [observable grid source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/grid/grid_obs.cc)

## function list:

### base class

- **``Grid::Grid``**
```
# input arguments
-
# return
-
```
> Default ``Grid`` constructor.

- **``Grid::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It serves as a virtual definition.

- **``Grid::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It serves as a virtual definition.

- **``Grid::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It serves as a virtual definition.

### regular magnetic grid

- **``Grid_breg::Grid_breg``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``read_permission`` or ``write_permission`` is true, it calls the ``Grid_breg::build_grid`` function.
The ``read_permission``, ``write_permission`` and ``build_permission`` are defined during parsing the external parameter file.

- **``Grid_breg::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It allocates the memory for regular magnetic field, according to the given sampling size.

- **``Grid_breg::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical regular magnetic field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_breg::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical regular magnetic field in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### random magnetic grid

- **``Grid_brnd::Grid_brnd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``build_permission`` or ``write_permission`` is true, it calls the ``Grid_brnd::build_grid`` function.
The ``build_permission`` is set as true whenever we need to make a realization of the field.

- **``Grid_brnd::~Grid_brnd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> This destructor is in charge of cleaning the FFTW plans and releasing memory consumption.

- **``Grid_brnd::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It not only allocates the memory required for conducting FFT but also prepares the plans for the transform.
In hammurabi X we use in-plan transform by default.

- **``Grid_brnd::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical random magnetic field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_brnd::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical random magnetic field in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### regular thermal electron grid

- **``Grid_tereg::Grid_tereg``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``read_permission`` or ``write_permission`` is true, it calls the ``Grid_tereg::build_grid`` function.

- **``Grid_tereg::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It allocates the memory for regular thermal electron field, according to the given sampling size.

- **``Grid_tereg::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical regular thermal electron field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_tereg::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical regular thermal electron field in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### random thermal electron grid

- **``Grid_ternd::Grid_ternd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``build_permission`` or ``write_permission`` is true, it calls the ``Grid_ternd::build_grid`` function.
The ``build_permission`` is set as true whenever we need to make a realization of the field.

- **``Grid_ternd::~Grid_ternd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> This destructor is in charge of cleaning the FFTW plans and releasing memory consumption.

- **``Grid_ternd::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It not only allocates the memory required for conducting FFT but also prepares the plans for the transform.
In hammurabi X we use in-plan transform by default.

- **``Grid_ternd::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical random thermal electron field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_ternd::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical random thermal electron field in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### regular dust grid

- **``Grid_dreg::Grid_dreg``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``read_permission`` or ``write_permission`` is true, it calls the ``Grid_tereg::build_grid`` function.

- **``Grid_dreg::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It allocates the memory for regular dust field, according to the given sampling size.

- **``Grid_dreg::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical regular dust field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_dreg::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical regular dust in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### random dust grid

- **``Grid_drnd::Grid_drnd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``build_permission`` or ``write_permission`` is true, it calls the ``Grid_ternd::build_grid`` function.
The ``build_permission`` is set as true whenever we need to make a realization of the field.

- **``Grid_drnd::~Grid_drnd``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> This destructor is in charge of cleaning the FFTW plans and releasing memory consumption.

- **``Grid_drnd::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It not only allocates the memory required for conducting FFT but also prepares the plans for the transform.
In hammurabi X we use in-plan transform by default.

- **``Grid_drnd::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical random dust field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_drnd::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical random dust in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### cosmic ray electron grid

- **``Grid_cre::Grid_cre``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> If ``read_permission`` or ``write_permission`` is true, it calls the ``Grid_tereg::build_grid`` function.

- **``Grid_cre::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It allocates the memory for cosmic ray electron field, according to the given sampling size.

- **``Grid_cre::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out the numerical cosmic ray electron field in a grid-piont-wise manner into a binary file with its name defined in the parameter-set.

- **``Grid_cre::import_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It reads in the numerical cosmic ray electron field in a grid-piont-wise manner from a binary file with its name defined in the parameter-set.

### observable grid

- **``Grid_obs::Grid_obs``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It simply calls the ``build_grid`` function.


- **``Grid_obs::build_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It allocates memory for observable maps upon requests in the HEALPix form.

- **``Grid_obs::export_grid``**
```
# input arguments
const Param * (parameter-set)
# return
-
```
> It writes out simulated maps, during which, the emission will be rescaled from cgs units into the default units, e.g., synchrotron/dust emission maps are prepared in $$mK$$ CMB temperature, Faraday depth maps are prepared in $$rad/m^{2}$$, and dispersion measure maps are prepared in $$pc/cm^3$$.

## method list:

- **grid data arrangement**

> By default we prepare/define hyper-rectangular grids for storing the descrete numerical fields.
For a 3D field, which means the field it self is defined with respect to the spatial position,
the corresponding grid is designed with the number of sample points and total physical distance coverage.
While for a 4D field, e.g., the CRE spatial-spectral distribution, the spectral coordinate is logarithmic
(the energy $$E_i$$ at coordinate index $$i$$ reads $$E_i = E_0 \exp{i \delta_E}$$).

> We use array of double for storing field components, e.g., a numerical scalar field is hosted by a single array
with dimensional (local) indices following the x-y-z direction ordering and so the storage (global) index is calculated by
$$i_\mathrm{global} = n_y n_z i_x + n_z i_y + i_z ~,$$
where the direction definition is only important for how we attach the field on the grid and the observer's position.
In the 4D field case, the global index is calculated by
$$i_\mathrm{global} = n_y n_z n_E i_x + n_z n_E i_y + n_E i_z + i_E ~.$$
Setting the $$E$$ coordinate as the last coordinate can save memory accessing time consumption.

- **3D/4D linear interpolation**

> The linear interpolation is actually dimension-free, by which we mean, in each dimension (or direction),
the interpolation is defined linearly by knowing the distance of the target position with respect to the lower
and upper limit of the elemental volume in which the target position resides.
It also doesn't matter which direction is interpolated first.
For the 3D case, we use the method explicitly called tri-linear interpolation, while in the 4D case, the method
is adjusted but still follows the same idea.
