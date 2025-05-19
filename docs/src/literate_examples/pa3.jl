# # pa3.jl

using Ferrite: Lagrange, RefTriangle
using ClosedWaveguideDispersion

# We define the parameters
function n(x)
    return 1.0 + 0.5 * cos(2π * x[1]) * sin(2π * x[2])
end;
p = 1.0;
h = 1.0;
N = 150;

## Set up the grid 
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