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
aas-doi: 10.3847/1538-4365/ab72a2
aas-journal: Astrophysical Journal Supplement Series
---

# Summary

Realistic models of Galactic emission are critical ingredients in the study of both the 
multi-phase interstellar medium (ISM) and precision cosmology. 
Models for the polarized synchrotron emission, absorption, and Faraday rotation are 
required in both cases. 
For ISM studies, these trace the physical conditions in the Galaxy and this emission is a 
significant foreground signal for cosmological surveys. 
The fundamental physical principles of the radiative transfer processes have been well 
understood for around half a century [@Rybicki:1979], but simple analytic models are not 
sufficient to capture the local structure and non-linear phenomena revealed by the 
precision and range of modern measurements.

To meet the growing need for numerical simulation of the Galactic emission, ``hammurabi`` 
[@Waelkens:2009] was developed to simulate Galactic observables based on a 3D 
modelling of the  physical components of the Galaxy.
The original code design was, however, not up to the coding and numerical standards 
required for studies of the ISM and cosmic microwave background (CMB) foregrounds.
Furthermore, the focus in the Galactic emission modelling has migrated recently 
from assuming simple regular field structure to more complicated fields with turbulence, 
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
