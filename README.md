# Rigorous integration in Parabolic PDEs on the torus


Table of contents:


* [Introduction](#introduction)
* [Equations under study](#equations-under-study)
   * [The 2D Navier-Stokes equations](#the-2D-navier-stokes-equations)
   * [The Swift-Hohenberg PDE](#the-swift-hohenberg-pde)
   * [The Kuramoto-Sivashinskig PDE](#the-kuramoto-sivashinski-pde)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2403.10450) : "Recent advances about the rigorous integration of parabolic PDEs via fully spectral Fourier-Chebyshev expansions", M. Cadiot and J-P. Lessard, [ArXiv Link](https://arxiv.org/abs/2403.10450).

It provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# Equations under study

We present some numerical details for the applications presented in Section 4 in [[1]](https://arxiv.org/abs/2403.10450). Each PDE is defined on the torus (in space) and we prove constructively the existence of solutions to IVPs for a given initial data. 

## The 2D Navier-Stokes equations
The 2D Navier-Stokes equation for the vorticity reads 

$$\omega_t = \nu \Delta \omega + \left(\Delta^{-1}\nabla \times \omega\right)\cdot \nabla w$$

and is associated to an initial condition $\omega(0,x) = b(x)$. In the code "NS_main.jl", we choose $b(x) = 2\sin(2x_1)\sin(x_2) - 2 \sin(x_1)\sin(2x_2) + 0.8\sin(x_1)\sin(3x_2) - 0.8 \sin(3x_1)\sin(x_2)  + 1.2\sin(2x_1)\sin(3x_2) - 1.2 \sin(3x_1)\sin(2x_2)$ and $\nu = 0.2$. We first take a time step of size $h=0.03$ and solve 73 consecutive IVPs of length 0.03 (starting at $b$). The propagation of error from one integration to the next is described in the Appendix of  [[1]](https://arxiv.org/abs/2403.10450). 

In this proof we choose $A = I_d - \mathcal{K}$, where $\mathcal{K}$ is the compact operator defined in Section 2 of  [[1]](https://arxiv.org/abs/2403.10450). In particular, such a simplication allows to compute simple bounds for Theorem 2.2 and we get
- $Y = (1+ ||\mathcal{K}||)||F(\overline{U})||$
- $Z1 = ||\mathcal{K}||^2$
- $Z2 = (1 + ||\mathcal{K}||) ||\mathcal{D}_{C_0} DQ(\overline{U})||$.

In particular, we have that $||\mathcal{K}|| = \mathcal{O}(h)$, so since the time step $h$ is small, we still obtain that $Z1<1$ and we are able to verify the hypotheses of Theorem 2.2.

## The Swift-Hohenberg PDE

The 1D Swift-Hohenberg equation reads 

$$ u_t = (\alpha-1)u - 2u_{xx} - u_{xxxx} - u^3$$

and is associated to an initial condition $u(0,x) = b(x)$. In the code "SH_main.jl", we choose $b(x) = 0.02 cos(x)$ and $\alpha = 8.1$. We integrate the PDE for a time step of size $h=3$. In particular, the approximate solution $U_0$ is stored in the file "U0_SH.jld2" and contains a sequence of Fourier-Chebyshev coefficients with an order $K0 = 24$ in Fourier and $N0 = 150$ in Chebyshev.

## The Kuramoto-Sivashinski PDE

The 1D Kuramoto-Sivashinski reads 

$$u_t = -u_{xx} - \alpha u_{xxxx} + \nu uu_x$$

and is associated to an initial condition $u(0,x) = b(x)$. In the code "KS_main.jl", we choose $b$ to be a point on a periodic orbit, $\nu = -2$ and $\alpha = 0.127$. We integrate the PDE for a time step of size $h= 2.244333563761175$, corresponding to the period of the periodic orbit. In fact, we demonstrate the ability to integrate for a time step being a period. The approximate solution $U_0$ is stored in the file "U0_KS.jld2" and contains a sequence of Fourier-Chebyshev coefficients with an order $K0 = 31$ in Fourier and $N0 = 80$ in Chebyshev.


 # Utilisation and References

The files C0_i.jld2 and C1_i.jld2 ($i \in \{1, \dots, 10\}$) contain the precomputed values of $C_0(\mu)$ and $C_1(\mu)$ described in Section 3.1 of [[1]](https://arxiv.org/abs/2403.10450). These files have to be in the current folder and will be used directly in the proof. In particular, such values have been computed thanks to the code "computation_constants.jl".

Each "main" code is associated to a list of functions, which are necessary for the rigorous computations of the proof. The 3 lists essentially contains the same functions, and only differ at some very specific points. Such differences arise from the optimization of the given problem. The code corresponding to the Kuramoto-Sivashinski PDE contains a more detailled description of the numerical techniques. 


 The code is build using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) 
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl)
 - [MATLAB](https://github.com/JuliaInterop/MATLAB.jl)
 
 
 # License and Citation
 
This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{IVP-Cheb.jl,
  author = {Matthieu Cadiot},
  title  = {{IVP-Cheb}.jl},
  url    = {https://github.com/matthieucadiot/IVP-Cheb.jl},
  note = {\url{ https://github.com/matthieucadiot/IVP-Cheb.jl},
  year   = {2025}
}
```


# Contact

You can contact me at :

matthieu.cadiot@mail.mcgill.ca
