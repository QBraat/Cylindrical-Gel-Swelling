# Cylindrical-Gel-Swelling
Code for spring network model and poroelastic models (M.Sc project Quirine Braat)

Code repository for the M.Sc. project of Quirine Braat.

Folders contain the different codes that have been used to generate the results for the graduation thesis 'Anisotropic Swelling of a Cylindrical Hydrogel'. The individual folders contain README files with specific information regarding their contents. 

The notebooks have been used to generate data, and figures as part of the graduation project. The corresponding data and further details are available on Zenodo with doi10.5281/zenodo.5495244. 

_________________________________________________________________________________________________

Some general comments about the code:


The codes are divided in three folders, namely
(1) 'Cylindrical Poroelastic Model' containing the code of the cylindrical poroelastic model and its corresponding analysis codes
(2) 'Spherical Poroelastic Model' containing the code of the spherical poroelastic model and its corresponding analysis codes
(3) 'Spring models' containing the code of the spring models and the corresponding analysis


For the poroelastic models, the data can be generated with Cylinder-finalversion-annotated.wl (cylinder) or Spherical-swelling-notebook.wl (sphere). The output can be analyzed with other notebooks provided in the corresponding folders.

If one wants to implement the transition, one must uncomment the indicated line to the module ParametersLoop[]. 

The analysis/comparison codes are based on the generated data during the project. With the actual data (which can be retrieved from ZENODO), it is possible to generate the plots yourself. One can also generate new data themselves. Here, the codes used for analysis are provided as an example of how the data can be compared and/or analyzed. These files can be used as a starting point for analysis of the generated data.

The original notebooks are built in Mathematica version 12.1. 
