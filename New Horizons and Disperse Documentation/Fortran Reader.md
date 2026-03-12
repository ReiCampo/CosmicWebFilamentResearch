
### _init_.py 
The role of the __init__.py function is to take Fortran code and translate it in a way that we can work with Fortran files in python. 

The function has several classes that can:
- Read in the Halo Shapes calculated from the outputs of the catolog files found on the infinity server.
- Read Treebricks information from HAGN and New Horizons. These classes read low precision Treebricks.
- Read in New Horizon's Treebricks from the W0 haloshape dictionary.
- Read in High Precision Treebricks from New Horizons. 
- Creates a galaxy dictionary out of the binary files from Fortran outputs on Infinity.
- Reads in data about the Fortran data cube that was used in the simulation. 
- Reads in filament information from a filament dictionary that came from the ASCII NDSKL file in 2D and 3D.

Each class will be documented in greater detail on different trees