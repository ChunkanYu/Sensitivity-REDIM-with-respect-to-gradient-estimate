# Sensitivity-REDIM-with-respect-to-gradient-estimate

The MatLAB code provided in this project is to deal with the sensivitiy analyses of Reaction-Diffusion Manifolds (REDIMs) with respect to the gradient estimate. 
A simple test example is investigated here:

<img src="https://render.githubusercontent.com/render/math?math=\begin{array}{rllll}
 \frac{\partial\psi_1}{\partial t} &
 =- k_1 \psi_1&
 +d \frac{\partial^2\psi_1}{\partial x^2}  &\qquad \psi_1(0) = 1      &\psi_1(1) = 0 
 \\
  \frac{\partial\psi_2}{\partial t} &
  =\phantom{-}k_1 \psi_1- k_2  \psi_2 &
  +d \frac{\partial^2\psi_2}{\partial x^2}  &\qquad \psi_2(0) = 0      &\psi_2(1) = 0 
 \\
  \frac{\partial\psi_2}{\partial t} &
  = \phantom{-}k_2  \psi_2 - k_3  \psi_3 &
  +d \frac{\partial^2\psi_3}{\partial x^2}  &\qquad \psi_3(0) = 0      &\psi_3(1) = 0
 \end{array}">


