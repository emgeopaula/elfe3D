# elfe3D with updated features to calculate sensitivities for inversion, to become v2.0.0
Modelling with the total **el**ectric field approach using **f**inite **e**lements in **3D**

!!! Under development for enabling inversion !!!

*Output for inversion:*

- model: find out in which format required

for testing use 1D real array of log10 transformed resistivities, has to be also read in eventually and backtransformed to resistivities:

inv_model = log10(10) = 1.0 (size 43 with q1.8 model)

- forward data as forward_data (1D real array) from subroutine order_forward_data in mod_calculate_tf.f90

for testing:

with current input filet (three receivers), only Ex component written to forward_data

forward_data = 4.44562188E-05   2.80089204E-05   1.97479819E-05   3.22833081E-09   2.76825052E-09   2.51720481E-09

- sensitivities as Jvec (1D real array), JTvec (1D real array), Jrows (1D integer array), Jcols (1D integer array) from subroutine compute_Jvec_JTvec in mod_sensitivities.f90

for testing:
   Jvec = 666.0_dp
   JTvec = 444.0_dp
   Jrows = 6
   Jcols = 4

_About:_

`elfe3D` is a 3D forward modelling code that can simulate electric and magnetic field responses from frequency-domain controlled-source electromagnetic geophysical setups. It uses tetrahedral meshes and first-order finite-element approximations. In addition, adaptive mesh refinement approaches are implemented.

_Statement of need:_

`elfe3D`  solves forward problems arising from the curl-curl equation in terms of the total electric field using a direct forward solver. The code is designed for Earth Scientists who want to simulate electric and magnetic field responses originating from a transmitter and the interaction of its transmitted signal with the 3D Earth. This so-called controlled-source electromagnetic method is used to search for resources and environmental applications, such as geothermal energy, minerals or groundwater. The air and the Earth’s subsurface consist of cells hosting variable model parameters: isotropic electric resistivities and magnetic permeabilities. Compared to standard electromagnetic geophysical simulation software, `elfe3D` excels in flexibility regarding subsurface geometries and survey settings, i.e. receivers can be arbitrarily placed in the modelling domain and the electrical properties can be flexibly distributed in the subsurface upon model design. Implemented adaptive mesh refinement approaches can automatically design problem-specific meshes and optimise computational load and solution accuracy.

_Contributions:_

An earlier version of the code that `elfe3D` is based on was developed by Paula Rulff with contributions from Laura Maria Buntin and Thomas Kalscheuer at Uppsala University from 2018-2023 financed by the Smart Exploration project (European Union’s Horizon 2020 funding, grant agreement No. 775971).

The present version of `elfe3D` was released in 2024 under the Apache License, Version 2.0. Further developments of `elfe3D` by Paula Rulff, now at Delft University of Technology, are ongoing. Suggestions for improvements are welcome!

If you would like to report bugs in `elfe3D`, suggest specific ideas for improvement or seek support, please open an issue or send an email to p.rulff@tudelft.nl.

If you would like to contribute to `elfe3D`, please open a pull request or send an email to p.rulff@tudelft.nl. Upon preparing your contribution, please check that the code compiles, run the provided example test and compare your results to the reference solutions. Update this `README.md` by including an overview of the changes that you made. Update the manual `elfe3D/elfe3D/README.md` by including a description of the new features you implemented and add specifications of new input parameters, if needed.

_Getting started:_

You find the `elfe3D` source code in `elfe3D/elfe3D/` and the manual including instalation instructions in `elfe3D/elfe3D/README.md`.
`elfe3D` can be compiled with the provided Makefile.
Note that, the open source mesh generator `tetgen` and the direct solver `MUMPS` must be installed additionally. 

Be aware that some `tetgen` versions are not working properly. The test example in `elfe3D/elfe3D/in` can be used to test if your tetgen version is working. Its should run without warnings, PLC errors or intersections. Check that your mesh is a closed 3D-cube with three regions (air, earth, anomaly) by visualising it, e.g. with `ParaView`.

We have encountered bugs when using most available version of tetgen. One way to correctly generate mesh files is to use the tetgen conda package v1.5.0. If you have conda or miniconda installed, run the following command to install it:

```
conda install -c conda-forge tetgen=1.5.0
```

_Tests:_

The mesh file of an example model is located in `elfe3D/elfe3D/in` and reference solutions for this example in `elfe3D/elfe3D/out`. You can use them to test, if the code runs properly and produces the expected results.


_Credits:_

If you publish results generated with `elfe3D`, please give credit to the `elfe3D` developers by citing:

Paula Rulff, Laura M Buntin, Thomas Kalscheuer, Efficient goal-oriented  mesh refinement in 3-D finite-element modelling adapted for controlled source electromagnetic surveys, Geophysical Journal International, Volume 227, Issue
3, December 2021, Pages 1624–1645, https://doi.org/10.1093/gji/ggab264

and refer to the `elfe3D` version you used via the ZENODO DOI: https://doi.org/10.5281/zenodo.13309721

Do not forget to acknowledge `MUMPS` and `TetGen` developers!

