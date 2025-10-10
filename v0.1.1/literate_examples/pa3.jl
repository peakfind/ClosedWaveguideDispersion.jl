# # PA-Three

# ## Problem

# PA-Three is an example in [ehrhardt-sun-zheng2009](@cite).
# + ``q(x_{1}, x_{2}) = 1 + 0.5 \cos(2\pi x_{1}) \sin(2\pi x_{2})``.
# + Neumann boundary condition on ``\partial \Omega``.

# ## Code

# Load the necessary packages
using Ferrite: Lagrange, RefTriangle
using ClosedWaveguideDispersion

# In this example, we only need to define a new refractive index.
function n(x)
    return 1.0 + 0.5 * cos(2π * x[1]) * sin(2π * x[2])
end;

# The parameters related to the periodic cell and the discrete Brillouin zone. 
p = 1.0;
h = 1.0;
N = 150;

# Set up the grid 
grid = setup_grid(lc=0.05, period=p, height=h)

## Define the interpolation 
ip = Lagrange{RefTriangle, 1}()

## Set up the FE values: CellValues 
cv = setup_fevs(ip)

## Set up the DofHandler 
dh = setup_dofs(grid, ip)

## Set the boundary conditions
cst = setup_bdcs(dh, period=p)

## Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

## Discretize the Brillouin zone
bz = collect(range(-π/p, 2π/p, N))

## Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=6)

## Plot the dispersion diagram
plot_diagram(bz, μ, period=p)

# ## Plain code
# ```julia
# using Ferrite: Lagrange, RefTriangle
# using ClosedWaveguideDispersion
# 
# function n(x)
#     return 1.0 + 0.5 * cos(2π * x[1]) * sin(2π * x[2])
# end
# 
# p = 1.0
# h = 1.0
# N = 150
# 
# grid = setup_grid(lc=0.05, period=p, height=h)
# 
# ip = Lagrange{RefTriangle, 1}()
# cv = setup_fevs(ip)
# dh = setup_dofs(grid, ip)
# 
# cst = setup_bdcs(dh, period=p)
# 
# A = allocate_matries(dh, cst)
# B = allocate_matries(dh, cst)
# 
# bz = collect(range(-π/p, 2π/p, N))
# μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=6)
# 
# plot_diagram(bz, μ, period=p)
# ```