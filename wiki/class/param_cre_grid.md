## member list of **``param_cre_grid``**

- **``filename``**
> Binary file name of numerical cosmic ray electron field input.

- **``build_permission``**
> Permission for building the grid for cosmic ray electron field.

- **``read_permission``**
> Permission for reading in external cosmic ray electron field.

- **``write_permission``**
> Permission for writing out internal cosmic ray electron field.

- **``x_max``**
> Upper limit of x direction grid physical range.

- **``x_min``**
> Lower limit of x direction grid physical range.

- **``y_max``**
> Upper limit of y direction grid physical range.

- **``y_min``**
> Lower limit of y direction grid physical range.

- **``z_max``**
> Upper limit of z direction grid physical range.

- **``z_min``**
> Lower limit of z direction grid physical range.

- **``E_max``**
> Upper limit of E (energy) direction grid physical range.

- **``E_min``**
> Lower limit of E (energy) direction grid physical range.

- **``E_fact``**
> E (energy) direction logarithmic difference, which can be uderstood by
$$E_\mathrm{max} = E_\mathrm{min} \exp\{n_E E_\mathrm{fact}\} ~,$$
where $$n_E$$ is given by ``nE``.

- **``nx``**
> Number of grid vertices (or samples) in x direction.

- **``ny``**
> Number of grid vertices (or samples) in y direction.

- **``nz``**
> Number of grid vertices (or samples) in z direction.

- **``nE``**
> Number of grid vertices (or samples) in E (energy) direction.

- **``cre_size``**
> Full number of grid vertices (or samples), equals to the multiplication of ``nx``, ``ny``, and ``nz``.
