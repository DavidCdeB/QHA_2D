#

Merging `Automating_QHA_plots` and `QHA_step_clean` 

## Table of Contents

- [Install](#install)



## What is the QHA program ?
 
 `QHA` is a program for computational chemistry and physics that performs the quasi-harmonic approximation reading the frequencies at each volume calculated with [CRYSTAL](http://www.crystal.unito.it/index.php). 
 
1. Extracts all the frequencies within all the **k** points in the supercell

[comment]: <> ( * Fits the frequency of each normal mode with respect to the volume.)

2. Calculates the pressure at finite temperature, as well as the entropy, Helmholtz and Gibbs free energy.
3. Produces tables summarizing the results for all the volumes analyzed.
4. For a pair of solid state phases, if produces a PDF document with plots of the following type: <a href="https://www.codecogs.com/eqnedit.php?latex=F(V;T),&space;\,&space;\,&space;P(V;T),&space;\,\,&space;G(P;T)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F(V;T),&space;\,&space;\,&space;P(V;T),&space;\,\,&space;G(P;T)" title="F(V;T), \, \, P(V;T), \,\, G(P;T)" /></a>:

5. Outputs the pressure-temperature phase diagram for the thermodynamic phase stability of both solid phases: 

<img src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/PT_phase_Boundary_edit.png" width="256" height="256" title="Github Logo">


 
The program was developed as part of [David Carrasco de Busturia PhD project](https://www.imperial.ac.uk/people/d.carrasco-de-busturia/) at [Prof. Nicholas Harrison's Computational Materials Science Group](http://www.imperial.ac.uk/computational-materials-science/), Imperial College London. The program was used to investigate the phase diagram and phase transitions mechanisms on the calcium carbonate system.
