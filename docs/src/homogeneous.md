```@meta
CurrentModule = ClosedWaveguideDispersion
```

# Homogeneous case with Neumann boundary condition

```julia
# case1.jl
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