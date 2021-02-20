<!-- [![Paper](https://img.shields.io/badge/paper-arXiv%3A1912.09947-B31B1B.svg)](https://arxiv.org/abs/1912.09947) -->
<!-- [![DOI](https://zenodo.org/badge/214220909.svg)](https://zenodo.org/badge/latestdoi/214220909) -->



# Effective Two-Dimensional Bose–Hubbard Model for Helium on Graphene
[Jiangyong Yu](https://github.com/yuj120254), Ethan Lauricella, Mohamed Elsayed, Kenneth Shepherd Jr., [Nathan S. Nichols](https://github.com/nscottnichols), Todd Lombardi, Sang Wook Kim, Carlos Wexler, [Juan M. Vanegas](http://vanegaslab.org/), Taras Lakoba, Valeri N. Kotov, and [Adrian Del Maestro](https://delmaestro.org/adrian)

<!-- [arXiv:1912.09947](https://arxiv.org/abs/1912.09947) -->

### Abstract
An exciting development in the field of correlated systems is the possibility of realizing two-dimensional (2D) phases of quantum matter. For a systems of bosons, an example of strong correlations manifesting themselves in  a 2D environment is provided by helium adsorbed on graphene.  We construct the effective Bose-Hubbard model for this system which involves hard-core bosons (U=∞), repulsive nearest-neighbor (V>0) and small attractive (V'&#60;0) next-nearest neighbor interactions. The mapping onto the Bose-Hubbard model is accomplished by a variety of many-body techniques which take into account the strong He-He correlations on the scale of the graphene lattice spacing. Unlike the case of dilute ultracold atoms, where interactions are effectively point-like, the detailed microsocpic form of the short range electrostatic and long range dispersion interactions in the helium-graphene system are crucial for the emergent Bose-Hubbard description.  The result places the ground state of the first layer of <sup>4</sup>He adsorbed on graphene deep in the commensurate solid phase with 1/3 of the sites on the dual triangular lattice occupied. Because the parameters of the effective Bose-Hubbard model are very sensitive to the exact lattice structure, this opens up an avenue to tune quantum phase transitions in this solid-state system.

### Description
This repository includes information, code, scripts, and data to generate the figures in a paper.

### Requirements
The [raw data](https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/tree/master/data) in this project was generated via the various methods described in the text.  This includes quantum Monte Carlo, density functional theory, etc.  To regenerate all quantum Monte Carlo data from scratch, information on obtaining the source code can be found [here](https://code.delmaestro.org).

* QMC Data: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4553524.svg)](https://doi.org/10.5281/zenodo.4553524)
* [graphenetools-py](https://github.com/nscottnichols/graphenetools-py)
* [pimcscripts](https://github.com/DelMaestroGroup/pimcscripts)
* [dgutils](https://github.com/DelMaestroGroup/dgutils)
* [heprops](https://github.com/agdelma/heprops)

### Support
This work was supported, in part, under NASA grant number 80NSSC19M0143.
Computational resources were provided by the NASA High-End Computing (HEC) Program through the NASA Advanced Supercomputing (NAS) Division at Ames Research Center. 

[<img width="100px" src="https://www.nasa.gov/sites/all/themes/custom/nasatwo/images/nasa-logo.svg">](https://www.nasa.gov/offices/education/programs/national/epscor/home/index.html)

### Figures

#### Figure 01: Adsorption Cell
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/graphene_cell.png" width="400px">

#### Figure 02: Dual Triangular Lattice
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/triangular.svg" width="400px">

#### Figure 03: Mean Field Phase Diagram
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/phase_diagram.svg" width="400px">

#### Figure 04: Interaction and Adsorption Potentials
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/He_Graphene_Potential.png" width="400px">

#### Figure 05: Adsorbed 4He: z-Wavefunctions
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/shooting_method_rho.svg" width="400px">

#### Figure 06: Adsorbed 4He: Band Structure
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/band_structure.svg" width="400px">

#### Figure 07: Adsorbed 4He: xy-Wavefunctions
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/wannier_cut.svg" width="400px">

#### Figure 08: Adsorbed 4He: Tunneling Potential
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/tunneling_steele.svg" width="400px">

#### Figure 09: Single Particle Properties
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/N_eq_1_density.svg" width="400px">

#### Figure 10: V'
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/density_n_1o3_NA_48.svg" width="400px">

#### Figure 11: Layering Analysis
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/rho_z_Lz.svg" width="400px">

#### Figure 12: Best fit analysis for confinement box size
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/chi2_vs_Lz.svg" width="400px">

#### Figure 13: Equation of State
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/EoS_Teq0.svg" width="400px">

#### Figure 14: V and V'
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/density_n_1_NA_48.svg" width="400px">

#### Figure 15: Extracting V from QMC
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/VBH_vs_Lz.svg" width="400px">

#### Figure 16: &beta; and &tau; scaling
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/beta_tau_scaling.svg" width="400px">

#### Figure 17: Finite Size Scaling of V'
<img src="https://github.com/DelMaestroGroup/papers-code-BoseHubbardModelHeAdsorptionGraphene/blob/master/plots/Vprime_QMC_fss.svg" width="400px">





