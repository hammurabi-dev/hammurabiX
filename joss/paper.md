---
title: 'hammurabi X: a C++ package for simulating Galactic emissions'
tags:
- C++
- astronomy
- Galactic emission
authors:
- name: Jiaxin Wang
  orcid: 0000-0002-7384-7152
  affiliation: "1, 2, 3"
- name: Tess R. Jaffe
  orcid: 0000-0003-2645-1339
  affiliation: "4, 5"
- name: Torsten A. Ensslin
  orcid: 0000-0001-5246-1624
  affiliation: "6"
- name: Joe Taylor
  affiliation: "4"
affiliations:
- name: Scuola Internazionale Superiore di Studi Avanzati, Via Bonomea 265, 34136 Trieste, Italy
  index: 1
- name: Department of Astronomy, Shanghai Jiao Tong University, 800 Dongchuan Road, 200240 Shanghai, China
  index: 2
- name: Istituto Nazionale di Fisica Nucleare, Sezione di Trieste, Via Bonomea 265, 34136 Trieste, Italy
  index: 3
- name: Department of Astronomy, University of Maryland, College Park, MD, 20742, USA
  index: 4
- name: CRESST, NASA Goddard Space Flight Center, Greenbelt, MD 20771, USA
  index: 5
- name: Max Planck Institute for Astrophysics, Karl-Schwarzschild-Str. 1, D-85741 Garching, Germany
  index: 6
date: 6 October 2019
bibliography: paper.bib
aas-doi: 10.3847/xxxxx
aas-journal: Astronomical Journal
---

# Summary

Understanding the Galactic emission is critical for studying not only the multi-phase 
interstellar medium (ISM), but also for detailed investigations of the cosmic microwave 
background (CMB) radiation or the 21cm cosmology.
Both areas recognize the importance of physical modellings of the Galactic polarized 
synchrotron emission, absorption, and Faraday rotation. 
For ISM studies, these trace the physical conditions in the Galaxy, while for CMB studies, 
these provide the most important CMB foreground.

The fundamental physical principles of the radiative transfer processes have been well 
understood for around half a century [@Rybicki:1979], however, with the growing precision 
and range of observations we are overwhelmed by various local structures and non-linear 
phenomena within the Galaxy.
This is slowing down conceptual and theoretical advancements in the mentioned research 
areas since the observables are no longer analytically calculable in a high-resolution and 
non-perturbative regime.

To meet the growing need of numerical simulation of the Galactic emission, ``hammurabi`` 
[@Waelkens:2009] was developed to simulate Galactic observables based on a 3D 
modelling of the  physical components of the Galaxy.
The original code design was, however, not matching modern coding and numerical 
standards required by the current scientific developments in ISM and CMB foreground 
modelling.
Furthermore, the focus in the Galactic emission modelling has migrated recently 
from assuming simple regular fields structure to more complicated fields with turbulence, 
enhancing the need for an accurate state-of-the-art simulation package for Galactic emission.

To extract information on the Galactic magnetic fields,
@Boulanger:2018 proposed a Bayesian method for parametric and non-parametric Galactic magnetic field inference. 
As Bayesian inference is computationally expensive, it has to rest on a high-performance 
simulation package. To meet these requirements, the hammurabi code has been upgraded 
to ``hammurabi X``, with a complete package redesign including modern coding 
standards and support for multi-threading.

Technically, ``hammurabi X`` performs an efficient line-of-sight radiative transfer integral 
through the simulated Galaxy model using a HEALPix-based nested grid to produce 
observables such as Faraday rotation measure and diffuse synchrotron 
emission in full Stokes I, Q and U, while taking into account beam and depth depolarization 
as well as Faraday effects.

The scientific aim of ``hammurabi X`` is to provide the ISM and CMB communities 
with a versatile numerical package to simulate Galactic emission.
The numerical framework can also be applied to other settings, for example, recently 
``hammurabi X`` provided extra-galactic Faraday rotation maps resulting from 
(reconstructed) primordial magnetic fields [@Hutschenreuter:2018].

# Acknowledgements

We acknowledge technical contributions from Theo Steininger,
inspiring discussions with Ellert van der Velden,
and the support from SISSA HPC service and the MHPC program.
We thank Andrea Zonca, Duncan Watts and Daniel Foreman-Mackey 
for reviewing the software with constructive suggestions.

# References
