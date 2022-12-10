# Sensitivity Analysis of the REDIM with respect to gradient estimate

The MatLAB code provided in this project is to deal with the sensivitiy analyses of Reaction-Diffusion Manifolds (REDIMs) with respect to the gradient estimate. 
A simple test example is investigated here:

```math
\frac{\partial\psi_1}{\partial t} =- k_1 \psi_1  + D \frac{\partial^2\psi_1}{\partial x^2} ~~~~~~~~~~~~~
```

```math
\frac{\partial\psi_2}{\partial t} = + k_1 \psi_1 - k_2 \psi_2  + D \frac{\partial^2\psi_2}{\partial x^2}
```

```math
\frac{\partial\psi_3}{\partial t} = + k_2 \psi_2 - k_3 \psi_3 + D \frac{\partial^2\psi_3}{\partial x^2}
```

The Codes are structured as follows:
* detailed_solution

    --> This folder consists of source code to calculate the detailed solution
    
* 1D_REDIM
    
    --> This folder consists of source code to generate 1D REDIM reduced chemistry and the sensitivity with respect to gradient estimate
    
* 2D_REDIM
    
    --> This folder consists of source code to generate 2D REDIM reduced chemistry and the sensitivity with respect to gradient estimate

More explanation can be found in the corresponding sub-folders

