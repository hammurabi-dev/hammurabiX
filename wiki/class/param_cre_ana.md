## member list of **``param_cre_ana``**

- **``alpha``**
> $$\alpha$$, the cosmic-ray electorn spectral index intrinsic value.

- **``beta``**
> $$\beta$$, the cosmic-ray electron spectral index cylindrical radius correction parameter.

- **``theta``**
> $$\theta$$, the cosmic-ray electron spectral index cylindrical height correction parameter.

- **``r0``**
> $$r_0$$, the cylindrical radial scaling factor.

- **``z0``**
> $$z_0$$, the height scaling factor.

- **``E0``**
> $$E_0$$, the cosmic-ray normalization energy scale.

- **``j0``**
> $$j_0$$, the cosmic-ray electron normalization flux at $$E_0$$.

The explicit ``analytic`` cosmic-ray electron flux distribution reads

$$\Phi(E,r,z) = j_0(\frac{\beta}{\beta_0})\gamma_0^{-idx_\odot}\gamma^{idx}\exp(\frac{r_\odot-r}{r_0})sech^2(\frac{z}{z_0})$$

$$\gamma_0 = E_0/m_e c^2$$

$$idx = -\alpha + \beta r + \theta z$$

$$idx_\odot = -\alpha + \beta r_\odot + \theta z_\odot$$

where E is the cosmic-ray electron energy in GeV, r is the cylindrical radius with respect to the galactic center, z is the cylindrical height with respect to the galactic disk.
