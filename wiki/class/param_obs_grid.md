## member list of **``param_obs_grid``**

- ``write_permission``
> It tells if any observable output is expected. 
Its default value is false.

- ``nside_dm``
> It defines the output map resolution of dispersion measure.

- ``nside_fd``
> It defines the output map resolution of Fraraday rotation.

- ``nside_sync``
> It is a vector of unsigned integers, since one simulation routine can provide synchrotron emission at various frequencies.
It defines the output map resolution of synchrotron emission.

- ``total_shell``
> It defines the number of total shperical shells in simulation.

- ``nside_shell``
> It defines the spherical resolution for each shell, from the inner shell to the outer shell.
It is a vector of unsigned integers.

- ``cut_shell``
> It defines the ratios of shell dividing radius with respect to the simulation radius limit ``oc_r_min`` and ``oc_r_max``.
It is a vector of doubles.

- ``radii_shell``
> It defines the physical shell dividing radii.
It is a vector of doubles.
It is calculated by using ``cut_shell``, ``oc_r_min`` and ``oc_r_max``.

- ``oc_r_min``
> It defines the minimal physical simulation spherical radius with respect to the observer.

- ``oc_r_max``
> It defines the maximal physical simulation spherical radius with respect to the observer.

- ``gc_r_min``
> It defines the minimal physical simulation spherical radius with respect to the galactic center.

- ``gc_r_max``
> It defines the maximal physical simulation spherical radius with respect to the galactic center.

- ``gc_z_min``
> It defines the lowest physical simulation z direction position with respect to the galactic center.

- ``gc_z_max``
> It defines the highest physical simulation z direction position with respect to the galactic center.

- ``oc_r_res``
> It defines the LoS integral resolution for all shells.

- ``do_dm``
> It tells if dispersion measure map is required as simulation output.
Its default value is false.

- ``do_fd``
> It tells if Faraday depth map is required as simulation output.
Its default value is false.

- ``do_sync``
> It tells if synchrotron emission map is required as simulation output.
Its default value is false.

- ``sim_fd_name``
> It defines the file name of Faraday depth output map.

- ``sim_dm_name``
> It defines the file name of dispersion measure output map.

- ``sim_sync_name``
> It defines the file name of synchrotron emission output map.

- ``do_mask``
> It tells if masking is required in observable outputs.
Its default value is false.

- ``mask_name``
> It defines the input mask map file name.
Note that the mask map is applied universally to all output observables.
