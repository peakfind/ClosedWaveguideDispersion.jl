```@meta
CurrentModule = ClosedWaveguideDispersion
```

# Homogeneous case with Neumann boundary condition

## Problem

This problem in given as an example in [zhang2021](@cite).
+ The function ``q(x_{1}, x_{2}) = 1`` is a constant function.
+ Neumann boundary condition on ``\partial \Omega``

## Code

First we should load ClosedWaveguideDispersion.jl for implemented functions. We also need Ferrite.jl to define the interpolation.
```example homo
using ClosedWaveguideDispersion
using Ferrite: Lagrange, RefTriangle
```
Since we consider the homogeneous waveguide, we need to define the refraction index.
```example homo
# Refractive index
function n(x)
   return 1.0 
end;

# Period of the waveguide
p = 1.0;
# Height of the waveguide
h = 1.0;
# Number of points in the discrete Brillouin zone
N = 100;
```
We use [`setup_grid`](@ref) to generate mesh for the periodic cell with period `p` and height `h`.
```example homo
# Set up the grid 
grid = setup_grid(lc=0.05, period=p, height=h)
```
Then we need to define the interpolation, `CellValues` and `DofHandler` which are needed by Ferrite.jl.
```example homo
# Define the interpolation 
ip = Lagrange{RefTriangle, 1}()

# Set up the FE values: CellValues 
cv = setup_fevs(ip)

# Set up the DofHandler 
dh = setup_dofs(grid, ip);
```
In this example, the Neumann boundary condition is satisfied naturally. So we only need to impose the periodic boundary condition in [`setup_bdcs`](@ref).
```example homo
# Set the boundary conditions
cst = setup_bdcs(dh, period=p)
```
Now we should prepare for the assembly of the eigenvalue problem.
```example homo
# Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

# Discretize the Brillouin zone
bz = collect(range(-π/p, π/p, N))
```
In [`calc_diagram`](@ref), `B` is generated only one time and `A` is generated with respect to ``\alpha`` in `bz`. Then we solve the generalized linear eigenvalue problem ``Ax = \mu Bx`` at each ``\alpha``.
```example homo
# Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=7)
```
Finally, we can plot our dispersion diagram.:
```example homo
# Plot the dispersion diagram
plot_diagram(bz, μ, period=p)
```

## Plain code
```julia
using ClosedWaveguideDispersion
using Ferrite: Lagrange, RefTriangle

# Parameters
function n(x)
   return 1.0 
end

p = 1.0
h = 1.0
N = 100

# Set up the grid 
grid = setup_grid(lc=0.05, period=p, height=h)

# Define the interpolation 
ip = Lagrange{RefTriangle, 1}()

# Set up the FE values: CellValues 
cv = setup_fevs(ip)

# Set up the DofHandler 
dh = setup_dofs(grid, ip)

# Set the boundary conditions
cst = setup_bdcs(dh, period=p)

# Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

# Discretize the Brillouin zone
bz = collect(range(-π/p, π/p, N))

# Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=7)

# Plot the dispersion diagram
plot_diagram(bz, μ, period=p)
```