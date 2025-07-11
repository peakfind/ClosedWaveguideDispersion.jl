# # Homogeneous case with Dirichelt boundary condition

# ## Problem
# 
# We consider homogeneous case with Dirichlet boundary condition on the boundaries of the closed waveguide.
# + The function ``q(x_{1}, x_{2}) = 1`` is a constant function.
# + homogeneous Dirichlet boundary condition on ``\partial \Omega``.

# ## Code 
# 
# We load the packages we need.
using Ferrite
using ClosedWaveguideDispersion

# Since we consider the homogeneous waveguide just like [Tutorial](@ref tutorial-homogeneous-neumann), we also define the refractive index as a constant.
function n(x)
    return 1.0 
end

# Similarly, here are the parameters related to the periodic cell and the discrete Brillouin zone.
p = 1.0;
h = 1.0;
N = 100;

# We need to implement a new function to impose the periodic boundary condition and the 
# Dirichlet boundary condition.
function my_bdcs(dh::DofHandler; period=1.0)
    cst = ConstraintHandler(dh)
    
    ## Periodic boundary condition on the left and right
    ## side of the periodic cell
    pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Ferrite.Vec{2}((period, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)
    
    ## Set Dirichlet boundary condition on the top and bottom boundaries of the periodic cell
    dfacets = union(getfacetset(dh.grid, "bottom"), getfacetset(dh.grid, "top"))
    dbc = Dirichlet(:u, dfacets, x -> 0)
    add!(cst, dbc)
    
    close!(cst)
    
    return cst 
end

# Set up the grid.
grid = setup_grid(lc=0.05, period=p, height=h)

# Basic steps needed by Ferrite.jl
ip = Lagrange{RefTriangle, 1}()
cv = setup_fevs(ip)
dh = setup_dofs(grid, ip)

## Set the boundary conditions
cst = my_bdcs(dh, period=p)

## Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

## Discretize the Brillouin zone
bz = collect(range(-π/p, π/p, N))

## Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=6)

## Plot the dispersion diagram
plot_diagram(bz, μ, period=p)

# ## Plain code 
# 
# ```julia
# @__CODE__
# ```