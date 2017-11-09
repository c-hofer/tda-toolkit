# tda-toolkit

This repository contains persistent homology related code. A final code release
is planned for late 2017.

# Installation

The `pershombox` package is dependent on some third party software tools which we do not provide here.
Each of those dependencies yields an executable which you have to copy to 
`pershombox/_software_backends/ext_lib`.

Where to find the executables/sources and how to name the executable is listed below. 
**Do not forget to `chmod` executables on Unix-based systems!**

1. `DIPHA`: [Source code](https://github.com/DIPHA/dipha). Rename to `dipha`.
2. `Perseus`: [Source code or precompiled executables](http://people.maths.ox.ac.uk/nanda/perseus/index.html). 
Rename to `perseus`.
3. `Dionysus`? (planned to support Dionysus 2)

# Main features
A short overview of the main features. For each of feature, there exists a tutorial in the 
`tutorials` subfolder.

### `toplex_persistence_diagrams`
Uses `Perseus` to calculate persistence diagrams of filtrated Toplex. [Tutorial](https://github.com/c-hofer/tda-toolkit/blob/tutorials_and_readme/tutorials/toplex_persistence_diagrams.ipynb)

### `cubical_complex_persistence_diagrams`
Uses `DIPHA` to calculate persistence diagrams of a filtrated cubical complex. [Tutorial](https://github.com/c-hofer/tda-toolkit/blob/tutorials_and_readme/tutorials/cubical_complex_persistence_diagrams.ipynb)

### `calculate_discrete_NPHT_2d`
Calculates a *normalized barycentric persistent homology transform* of a given binary 2D cubical complex.[Tutorial](https://github.com/c-hofer/tda-toolkit/blob/tutorials_and_readme/tutorials/discrete_2d_npht.ipynb)

### `calculate_discrete_NPHT_3d_Lebedev26`
Calculates a *normalized barycentric persistent homology transform* (residing on 
the 26-points Lebedev grid) of a given binary 3D cubical complex.
See [1].

### `distance_npht2D`
Calculates the 'shape' distance between two 2D persistent homology transforms. 

### `distance_npht3D_lebedev_26` 
Calculates the 'shape' distance between two 3D persistent homology transforms, proposed
in [1].

# References 
[[1]](http://wwwx.cs.unc.edu/~mn/sites/default/files/hofer2017_ipmi.pdf) 
C. Hofer, R. Kwitt, M. Niethammer, Y. Hoeller, E. Trinka and A. Uhl.    
*Constructing Shape Spaces from a Topological Perspective*, In: IPMI, 2017
```bash
@inproceedings{Hofer17a,
  author    = {C.~Hofer, R.~Kwitt, M.~Niethammer, Y.~Hoeller, E.~Trinka and A.~Uhl},
  title     = {Constructing Shape Spaces from a Topological Perspective},
  booktitle = {IPMI},
  year      = {2017}}
```


