#

<!--Merging `Automating_QHA_plots` and `QHA_step_clean`-->

## Table of Contents

- [What is the QHA program ?](#What is the QHA program ?)
- [What is the quasi-harmonic approximation ?](#What is the quasi-harmonic approximation ?)


## What is the `QHA_2D` program ?
 
 `QHA_2D` is a program for computational chemistry and physics that performs the quasi-harmonic approximation reading the frequencies at each volume calculated with [CRYSTAL](http://www.crystal.unito.it/index.php). 
 
* Extracts all the frequencies within all the **k** points in the supercell for a given volume.

<!--( * Fits the frequency of each normal mode with respect to the volume.)-->

* Calculates the pressure at finite temperature as a function of volume, as well as the entropy, Helmholtz and Gibbs free energy.
* Produces tables summarizing the results for all the volumes analyzed.
* For a pair of solid state phases, if produces a PDF document with plots of the following type: <a href="https://www.codecogs.com/eqnedit.php?latex=F(P;T),&space;\,&space;\,P(V;T),&space;\,&space;\,G(P;T),&space;\,&space;\,S(P;T)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F(P;T),&space;\,&space;\,P(V;T),&space;\,&space;\,G(P;T),&space;\,&space;\,S(P;T)" title="F(P;T), \, \,P(V;T), \, \,G(P;T), \, \,S(P;T)" /></a>

* Outputs the pressure-temperature phase diagram for the thermodynamic phase stability of both solid phases: 

<!--<img  align="center" src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/PT_phase_Boundary_edit.png" width="256" height="256" title="Github Logo"> -->
<p align="center">
  <img width="256" height="256" src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/PT_phase_Boundary_edit.png">
</p>

 
The program was developed as part of [David Carrasco de Busturia PhD project](https://www.imperial.ac.uk/people/d.carrasco-de-busturia/) at [Prof. Nicholas Harrison's Computational Materials Science Group](http://www.imperial.ac.uk/computational-materials-science/), Imperial College London. The program was used to investigate the phase diagram and phase transitions mechanisms on the calcium carbonate system.


## What is the quasi-harmonic approximation ?

Since the birth of quantum chemistry, almost every calculation was performed in the athermal limit (0K) and no pressure effects were considered (0Pa).
One of the most exciting challenges in an _ab intio_ calculation is to obtain information of the system at a finite temperature and pressure. This allow us to obtain a more realistic picture of the system in the everyday world, where temperature and pressure cannot be neglected and are indeed the driving force for many transformations in nature.

One of the most famous techniques for taking into account the effect of the temperature in the computed properties of molecules and crystals is _ab intio_ molecular dynamics [1-3], in which the Schrodinger equation is solved at each MD time step within the Born-Oppenheimer approximation. Unfortunetely, this is a very computationally expensive technique. Therefore, there is a huge interest in developing accurate and reliable models for the inclusion of a combined effect of pressure and temperature in the standard first-principles quantum chemical methods.

The sole effect of pressure can be easily taken into account
if the structure is fully relaxed through a volume constraint
geometry optimization process.
If this process is repeated over a series of equidistant volumes ranging
from compression to expansion
a energy-volume or pressure-volume relation equation of state is
obtained, which describes the behaviour of
a solid under compression and expansion in the athermal limit.
On the contrary, the inclusion of temperature is not so straightforward to implement.
Following the harmonic approximation (HA) formalism,
the lattice dynamics of crystal vibrations (i.e. phonons) are calculated at
each volume.

The HA indeed presents serious limitations that arise
from the fact that
only constant volume quantities
can be calculated, so that
thermal expansion cannot be explained.
Simpler and less computationally expensive,
but still effective tecnhiques can be used based on the so called quasi-harmonic
approximation (QHA), which
corrects the HA deficiencies
by maintaining the same harmonic expression
but introducing an explicit dependence of vibration phonon frequencies
on volume.
For a more detailed explanation, please check the pdf in this repository.
