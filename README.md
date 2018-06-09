#

<!--Merging `Automating_QHA_plots` and `QHA_step_clean`-->

# Table of Contents

<!-- - [What is the QHA program ?](#WhatisQHA)
- [What is the quasi-harmonic approximation ?](#Whatisquasi) -->
1. [What is the QHA program ?](#example)
2. [What is the quasi-harmonic approximation ?](#example2)
3. [Power of the quasi-harmonic approximation](#example3)
4. [Why is `QHA_2D` useful ?](#example4)
5. [Files needed for running `QHA_2D`](#example5)
6. [How to run `QHA_2D`](#example6)
7. [Test](#example7)
8. [How to cite](#example8)
9. [Contributing](#example9)
10. [References](#example10)


<a name="example"></a>
## What is the QHA program ?
 
 `QHA_2D` is a program for computational chemistry and physics that performs the quasi-harmonic approximation reading the frequencies at each volume calculated with [CRYSTAL](http://www.crystal.unito.it/index.php). 
 
* Extracts all the frequencies within all the **k** points in the supercell for a given volume.

<!--( * Fits the frequency of each normal mode with respect to the volume.)-->

* Calculates the pressure at finite temperature as a function of volume, as well as the entropy, Helmholtz and Gibbs free energy.
* Produces tables summarizing the results for all the volumes analyzed.
* For a pair of solid state phases, if produces a PDF document with plots of the following type: <a href="https://www.codecogs.com/eqnedit.php?latex=F(P;T),&space;\,&space;\,P(V;T),&space;\,&space;\,G(P;T),&space;\,&space;\,S(P;T)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F(P;T),&space;\,&space;\,P(V;T),&space;\,&space;\,G(P;T),&space;\,&space;\,S(P;T)" title="F(P;T), \, \,P(V;T), \, \,G(P;T), \, \,S(P;T)" /></a>

* Outputs the pressure-temperature phase diagram for the thermodynamic phase stability of both solid phases: 

<!--<img  align="center" src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/PT_phase_Boundary_edit.png" width="256" height="256" title="Github Logo"> -->
<p align="center">
  <img width="256" height="256" src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/PT_phase_Boundary.svg">
</p>

* The underlying criteria for producing this phase boundary is
by evaluating where two Gibbs free energy intersect:
<a href="https://www.codecogs.com/eqnedit.php?latex=G^{I}(P;T)&space;=&space;G^{II}(P;&space;T)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G^{I}(P;T)&space;=&space;G^{II}(P;&space;T)" title="G^{I}(P;T) = G^{II}(P; T)" /></a>

The program was developed as part of [David Carrasco de Busturia PhD project](https://www.imperial.ac.uk/people/d.carrasco-de-busturia/) at [Prof. Nicholas Harrison's Computational Materials Science Group](http://www.imperial.ac.uk/computational-materials-science/), Imperial College London. The program was used to investigate the phase diagram and phase transitions mechanisms on the calcium carbonate system.

<a name="example2"></a>
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

<a name="example3"></a>
## Power of the quasi-harmonic approximation

By implementing a quasi-harmonic approximation framework we are able to calculate the **phase diagram** of a solid substance (like the one shown in the above figure) and the **thermodynamic stability** of different polymorphs, without the need of running computationally expensive __ab initio__ molecular dynamics calculations. In addition, this leads to the exploration of **phase transitions** at a fnite temperature and pressure.


<a name="example4"></a>
## Why is `QHA_2D` useful ?

The actual version of CRYSTAL17 does perform an atomated quasi-harmonic approximation calculation for a given set of volumes. However, the QHA built-in module in CRYSTAL presents some deficiencies:

* If the phase transition is driven by a soft phonon mode, i.e. a phonon mode that becomes negative at a certain volume constraint, the `QHA` module in CRYSTAL will give some problems.

* Unfortunalety, the optimization at a constant volume is performed within the supercell scheme. If the supercell is big, (as it should be in order to ensure convergence the entropy), this might leads to an optimization porcess which is without doubt, more difficult in the supercell scheme: the cell is bigger, there are more atoms, and this can lead to convergence problems, or flase minima.

* Not relevant (unwanted) supercell phase transitions as a consequence of the optimization being performed in the supercell.
   
* Ideally, CRYSTAL should perform the optimization in the the primitive cell prior to making the supercell for the phonons calculation at a finite **k** point, but this is not so trivial to implement in the main code, according to the developers. Hopefully, this will be taken into account in future versions of the code. But for the moment, the `QHA_2D` code presented in this repository is an easy and effective solution for evaluating thermodynamic properties of crystals at a finite temperature and pressure (the real world).

<a name="example5"></a>
## Files needed for running `QHA_2D`
 

* Say you want to compute the pressure-temperature phase diagram of two
sold phases I and II.
 `QHA_2D` requires the frequency calculation outputs at each volume, for each of the two phases.
These frequencies calculations can be either in the Gamma point or at finite **k** points.

* The name of all these frequency outputs have to end as `*.out`

* Please ensure that you are using a sufficient big supercell for the entropy to be converged with the number of **k** points (mention the other code to sort out the supercell expansion matrix).

<a name="example6"></a>
## How to run `QHA`

* Get the code: `git clone https://github.com/DavidCdeB/QHA_2D`
* Give permissions to all the scripts: `chmod u+x *.sh *.py`
* Create the `Files_Outputs` folder inside the `QHA_2D` folder that has just been cloned: `cd ./QHA_2D && mkdir Files_Outputs`
* Create the folders that will contain the constant-volume frequency outputs for each phase: `mkdir Calcite_I && mkdir Calcite II`
* Copy all the frequencies outputs for each volume, for each phase, to the folders `Calcite_I` and `Calcite_II`. For example, `Calcite_I` folder will contain the frequency output for each `j`-th volume for the Calcite I phase.
* Remember that name of all these frequency outputs have to end as `*.out`
* The file system at this point looks like the following:

<p align="left">
  <img width="256" height="256" src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/file_system.svg">
</p>

The following is a summarized flow chart on the structure of codes:

<p align="center">
  <img src="https://github.com/DavidCdeB/QHA_2D/blob/master/Images_for_README_md/codes.svg">
</p>

* Run `./boundary_1_node.sh`

**_Prerequisites_**

To run, `QHA` requires Python with certain packages:

* Python 2.7 or higher.
    Packages: `numpy`, `scipy`, `re`, `os`, `glob`, `itertools`, `subprocess`, `sys` (All of these come with a default [Anacaonda](https://www.continuum.io/downloads) installation).

* Standard `bash` version in your system.

<a name="example7"></a>
## Test

Under the `TEST` folder, you will find all the programs
needed, together with a `Files_Outputs` folder with the frequency outputs of two phases: calcite I and calcite II.
If you run the program, you will obtain the `main.pdf` with all the plots needed.

<a name="example8"></a>
## How to cite

Please cite the following reference when using this code:

Carrasco-Busturia, D., Erba, A., Mallia, G., Mellan, T. A. and Harrison, N. M. "Computed phase stability and phase transition mechanisms in CaCO3 at finite temperature and pressure" _In progress_

<a name="example9"></a>
## Contributing

`QHA` is free software released under the Gnu Public Licence version 3. 
All contributions to improve this code are more than welcome.

* Have a look at GitHub's ["How to contribute"](https://guides.github.com/activities/contributing-to-open-source/#contributing).

* If you are familiar with `git`: fork this repository and submit a pull request.

* If you are not familiar with `git`: 

    * If something should be improved, open an issue here on GitHub
    * If you think a new feature would be interesting, open an issue
    * If you need a particular feature for your project contact me directly.

<a name="example10"></a>  
## References
  
  [1] R. Car, M. Parrinello, Phys. Rev. Lett. 1985, 55, 2471
  
  [2] F. Buda, R. Car, M. Parrinello, Phys. Rev. B 1990, 41, 1680
  
  [3] F. D. Vila, V. E. Lindahl, J. J. Rehr, Phys. Rev. B 2012, 85, 024303

