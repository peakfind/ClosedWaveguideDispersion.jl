```@meta
CurrentModule = ClosedWaveguideDispersion
```

# Tutorial

!!! note
    For simplicity, we restrict our discussion to the computation of dispersion diagrams. For more theoretical details, we refer to [photonic](@cite). 

## Problem

## Variational formulation

In the end, we conclude the important steps in the computation
+ We use the Finite element method to discretize the variational formulation. In ClosedWaveguideDispersion.jl, all Finite element codes are implemented by [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl)
+ After the Finite element discretization, we obtain a generalized linear eigenvalue problem parametered by ``\alpah``. We utilize [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl) to solve the generalized linear eigenvalue problems with fixed ``\alpha``.