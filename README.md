#

Merging `Automating_QHA_plots` and `QHA_step_clean` 

## Table of Contents

- [Install](#install)



## What is the QHA program ?
 
 `QHA` is a program for computational chemistry and physics that performs the quasi-harmonic approximation reading the frequencies at each volume calculated with [CRYSTAL](http://www.crystal.unito.it/index.php). 
 
 * Extracts all the frequencies within all the **k** points in the supercell
[comment]: <> ( * Fits the frequency of each normal mode with respect to the volume.)
 * Calculates the pressure at finite temperature, as well as the entropy, and Gibbs free energy
 * Produces tables summarizing the results for all the volumes analyzed.
 * Produces a PDF document with plots of the following type: 

<a href="https://www.codecogs.com/eqnedit.php?latex=F" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F" title="F" /></a>
 
The program was developed as part of [David Carrasco de Busturia PhD project] (http://www.imperial.ac.uk/computational-materials-science/people/) at [Prof. Nicholas Harrison's Computational Materials Science Group](http://www.imperial.ac.uk/computational-materials-science/), Imperial College London. The program was used to investigate the phase diagram and phase transitions mechanisms on the calcium carbonate system.
