# tda-toolkit

This repository contains persistent homology related code which can be used 
to implement the approaches from [1] and [2] (see [References](#references))

*A final code release is planned for late 2017.*

# Installation

The `pershombox` package is dependent on some third party software tools which we do not provide here.
In order to install `pershombox` you have to get those executables and tell `pershombox` where 
 to find them by editing the corresponding entry in 
 `pershombox/_software_backends/software_backends.cfg`.

Where to find the executables/sources and which `software_backends.cfg` entry corresponds to them 
is listed below. 
**Do not forget to `chmod` executables on Unix-based systems!**

1. `DIPHA`: [Source code](https://github.com/DIPHA/dipha). 
Entry: `dipha`.

2. `Perseus`: [Source code or precompiled executables](http://people.maths.ox.ac.uk/nanda/perseus/index.html). 
Entry: `perseus`.

3. `hera`: [Source code](https://bitbucket.org/grey_narn/hera.git). 
We need the `wasserstein_dist` executable in `geom_matching/wasserstein`. Entry: `hera_wasserstein_dist`.

We plan to also support [Dionysus (v2)](http://mrzv.org/software/dionysus2/) in the future.

## Exemplary DIPHA installation

```bash
git clone https://github.com/DIPHA/dipha.git
cd dipha
mkdir build
cmake ..
make -j4
```

Then manipulate the `software_backends.cfg`: 

```bash
[paths]
# Configure the paths to the backend software here
# e.g., dipha=/home/myHome/dipha
# do not forget to chmod +x on unix bases systems

dipha=<path/to/your/dipha/executable/here> # <-- This is your modification

hera_wasserstein_dist=

perseus=
```

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
[Tutorial](https://github.com/c-hofer/tda-toolkit/blob/tutorials_and_readme/tutorials/discrete_2d_npht.ipynb)

### `distance_npht3D_lebedev_26` 
Calculates the 'shape' distance between two 3D persistent homology transforms, proposed
in [1].

# References 
[[1]](http://wwwx.cs.unc.edu/~mn/sites/default/files/hofer2017_ipmi.pdf) 
C. Hofer, R. Kwitt, M. Niethammer, Y. Hoeller, E. Trinka and A. Uhl.    
**Constructing Shape Spaces from a Topological Perspective**, In: IPMI, 2017
```bash
@inproceedings{Hofer17a,
  author    = {C.~Hofer, R.~Kwitt, M.~Niethammer, Y.~Hoeller, E.~Trinka and A.~Uhl},
  title     = {Constructing Shape Spaces from a Topological Perspective},
  booktitle = {IPMI},
  year      = {2017}}
```

[[2]](https://arxiv.org/abs/1707.04041) 
C. Hofer, R. Kwitt, M. Niethammer and A. Uhl.     
**Deep Learning with Topological Signatures**, In: NIPS, 2017
```bash
@inproceedings{Hofer17b,
  author    = {C.~Hofer, R.~Kwitt, M.~Niethammer, and A.~Uhl},
  title     = {Deep Learning with Topological Signatures},
  booktitle = {NIPS},
  year      = {2017}}
```


