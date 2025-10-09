module ClosedWaveguideDispersion

using Gmsh
using Ferrite
using FerriteGmsh: togrid
using SparseArrays
using Arpack: eigs
using CairoMakie
using StaticArrays

include("fem.jl")
# For 1D periodicity
export setup_grid, setup_fevs, setup_dofs, setup_bdcs
export allocate_matries, assemble_A, assemble_B

# For 2D periodicity
export setup_grid_squareLattice
export assemble_A_TE, assemble_A_TM, assemble_B_TE, assemble_B_TM

include("evp_solver.jl")
# For 1D periodicity
export calc_diagram

# For 2D periodicity
export calc_diagram_TE, calc_diagram_TM

include("plot.jl")
export plot_diagram

include("lattice.jl")
export SquareLattice, IrreducibleBrillouin, get_discrete_irreducibleBrillouin

end
