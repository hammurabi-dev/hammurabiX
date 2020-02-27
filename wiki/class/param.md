## hammurabi X Param class

The ``Param`` class is designed for storing, parsing and using simulation related parameters in a collective way.
A parameter may be used in controlling the simulation resolution, output type, etc., or in defining galactic components.
Instead of presenting a method list, we provide a parameter list which would help users/developers/contributers to have a clear view.

- [header file](https://github.com/hammurabi-dev/hammurabiX/tree/master/include/param.h)
- [source file](https://github.com/hammurabi-dev/hammurabiX/tree/master/source/param/param.cc)

## function list:

- **``Param::Param``**
```
# input arguments
const std::string (parameter-set file name)
# output
-
```
> It parses the parameters from the XML type parameter file.

- **``Param::obs_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters (directly) related to observable outputs.

- **``Param::breg_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters related to regular magnetic fields.

- **``Param::brnd_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters related to random magnetic fields.

- **``Param::tereg_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters related to regular thermal electron fields.

- **``Param::ternd_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters related to random thermal electron fields.

- **``Param::cre_param``**
```
# input arguments
tinyxml2::XMLDocument * (parameter-set)
# output
-
```
> It initializes the parameters related to cosmic ray electron fields.

## parameter list:

### model independent parameters:

- [param_obs_grid](./param_obs_grid.md)
> Observable grid related parameters.

- [param_breg_grid](./param_breg_grid.md)
> Regular magnetic grid related parameters.

- [param_brnd_grid](./param_brnd_grid.md)
> Random magnetic grid related parameters.

- [param_tereg_grid](./param_tereg_grid.md)
> Regular thermal electron grid related parameter.

- [param_ternd_grid](./param_ternd_grid.md)
> Random thermal electron grid related parameters.

- [param_cre_grid](./param_cre_grid.md)
> Cosmic ray electron grid related parameters.

### (non-trivial) model dependent parameters:

- [param_breg_lsa](./param_breg_lsa.md)
> LSA (logarithmic-spiral-arm) magnetic field model parameters.

- [param_breg_jaffe](./param_breg_jaffe.md)
> Jaffe magnetic field model parameters.

- [param_brnd_global_es](./param_brnd_global_es.md)
> Ensslin-Steininger global random magnetic field model parameters.

- [param_brnd_local_mhd](./param_brnd_local_mhd.md)
> Parameterized MHD magnetic field model parameters.

- [param_tereg_ymw16](./param_tereg_ymw16.md)
> YMW16 thermal electron field model parameters.

- [param_cre_ana](./param_cre_ana.md)
> Analytic position-dependent-constant spectral index cosmic ray electron field model parameters.

### (testing) model dependent parameters:

- [param_breg_unif](./param_breg_unif.md)
> Uniform magnetic field model parameters.

- [param_tereg_unif](./param_tereg_unif.md)
> Uniform thermal electron field model parameters.

- [param_cre_unif](./param_cre_unif.md)
> Uniform cosmic-ray electron field model parameters.
