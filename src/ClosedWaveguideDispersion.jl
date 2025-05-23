module ClosedWaveguideDispersion

using Gmsh
using Ferrite
using FerriteGmsh: togrid
using SparseArrays
using Arpack: eigs
using CairoMakie

include("fem.jl")
export setup_grid, setup_fevs, setup_dofs, setup_bdcs
export allocate_matries, assemble_A, assemble_B

include("evp_solver.jl")
export calc_diagram

include("plot.jl")
export plot_diagram

end