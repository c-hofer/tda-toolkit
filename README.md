# tda-toolkit
My persistent homology related code. 
At the moment code release is planned until end of august. :) 


# Installation

The `pershombox` package is dependent on some third party software tools which I do not provide here.
Each of those dependencies yields an executable which you have to copy to 
`pershombox/_software_backends/ext_lib`.
Where to find the executables/sources and how to name the executable is listed below. 

1. `DIPHA`: [Source code](https://github.com/DIPHA/dipha). Name to `dipha`.
2. `Perseus`: [Source code or precompiled executables](http://people.maths.ox.ac.uk/nanda/perseus/index.html). 
Name to `perseus`.
3. `Dionysus`?

# Main features
A short overview of the main features. For each of them exists a tutorial in the `tutorials` subfolder.

###`toplex_persistence_diagrams`<b>
Uses `Perseus` to calculate persistence diagrams of filtrated Toplex. 

###`cubical_complex_persistence_diagrams`<b>
Uses `DIPHA` to calculate persistence diagrams of a filtrated cubical complex.

###`calculate_discrete_NPHT_2d`<b>
Calculates a *normalized barycentric persistent homology transform* of a given binary 2D cubical complex.

###`calculate_discrete_NPHT_3d_Lebedev26`<b>
Calculates a *normalized barycentric persistent homology transform* (residing on 
the 26-points Lebedev grid) of a given binary 3D cubical complex.
See [1].

###`distance_npht2D`
Calculates the 'shape' distance between two 2D persistent homology transforms. 

###`distance_npht3D_lebedev_26`
Calculates the 'shape' distance between two 3D persistent homology transforms. See [1].

# References 
[[1]](http://wwwx.cs.unc.edu/~mn/sites/default/files/hofer2017_ipmi.pdf)
Hofer C., Kwitt R, Niethammer M., Hoeller Y., Trinka E., Uhl A.: 
*Constructing Shape Spaces from a Topological Perspective*, 
IPMI, 
2017




