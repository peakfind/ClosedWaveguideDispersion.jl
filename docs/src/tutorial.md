```@meta
CurrentModule = ClosedWaveguideDispersion
```

# Tutorial

In this tutorial, we will present some basic knowledge of dispersion digrams and use an example to show how ClosedWaveguideDispersion.jl works.

## Problem

Let ``\Omega = \mathbb{R} \times (0, 1) \subset \mathbb{R}^{2}`` be a closed waveguide with the boundary
```math
\partial \Omega = \{x = (x_{1}, x_{2}) \in \mathbb{R}^{2} : x_{2} = 0 \text{ or } x_{2} = 1\}.
```
Suppose ``\Omega`` is filled by periodic medium with period ``p``. The medium can be characterized by the real-valued refactive index ``q(x_{1}, x_{2})`` satisfying
```math
q(x_{1} + p, x_{2}) = q(x_{1}, x_{2}), q \geqslant c > 0, \forall (x_{1}, x_{2}) \in \Omega.
```
Moreover, we assume that ``q(x_{1}, x_{2}) \in L^{\infty}(\Omega)``.

The wave propagation in the periodic closed waveguide is modeled by the following boundary value problem
```math
\begin{align*}
    &\Delta u + k^{2}q(x_{1}, x_{2})u = f \text{ in } \Omega, \\
    &\frac{\partial u}{\partial x_{2}} = 0 \text{ on } \partial \Omega,
\end{align*}
```
where ``f`` is a function in ``L^{2}(\Omega)`` with compact support, and ``k > 0`` is the real wavenumber.

## Floquet-Bloch theory

For simplicity, we only present necessary notions for computation. For more theoretical details, we refer to [fliss-joly2016](@cite) and [photonic](@cite).

The parameter ``\alpha`` is introduced by the Floquet-Bloch transform. And we always consider ``\alpha`` in the Brillouin zone ``(-\pi/p, \pi/p)``.

Transfering the dependence on ``\alpha`` from the function space to the PDE, we have the following boundary value problem
```math
\begin{align*}
&\Delta v + 2i\alpha \partial_{1}v + (k^{2}q(x_{1}, x_{2}) - \alpha^{2}) v = 0 \text{ in } \Omega_{0}, \\
&\frac{\partial v}{\partial x_{2}} = 0 \text{ on } \partial \Omega_{0}, \\
&v \text{ is periodic with respect to } x_{1}.
\end{align*}
```

## Variational formulation

In this section, we present the variational formulation of the boundary value problem: Find ``v \in H_{per}^{1}(\Omega_{0})`` satisfying
```math
\int_{\Omega_{0}} \nabla v \cdot \nabla \bar{\phi} - 2i\alpha \partial_{1} v \bar{\phi} - (k^{2}q(x_{1}, x_{2}) - \alpha^{2}) v \bar{\phi} dx = 0.
```
After finite element discretization, we can obtain a generalized linear eigenvalue problem
```math
\mathbf{A}_{\alpha} \mathbf{v} = k^{2} \mathbf{B} \mathbf{v},
```
where ``\mathbf{A}_{\alpha}`` comes from
```math
\int_{\Omega_{0}} \nabla v \cdot \nabla \bar{\phi} - 2i\alpha \partial_{1} v \bar{\phi} + \alpha^{2} v \bar{\phi} dx 
```
and ``\mathbf{B}`` comes from
```math
\int_{\Omega_{0}} q(x_{1}, x_{2}) v \bar{\phi} dx.
```

!!! note "important steps in the computation"
    + We use the Finite element method to discretize the variational formulation. In ClosedWaveguideDispersion.jl, all Finite element codes are implemented by [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl)
    + After the Finite element discretization, we obtain a generalized linear eigenvalue problem parametered by ``\alpha``. We utilize [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl) to solve the generalized linear eigenvalue problems with fixed ``\alpha``.