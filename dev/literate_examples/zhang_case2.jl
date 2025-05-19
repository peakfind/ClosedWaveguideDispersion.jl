# # Inhomogeneous case

using Ferrite: Lagrange, RefTriangle
using ClosedWaveguideDispersion

# Parameters
function n(x)
    r = sqrt(x[1]^2 + (x[2] - 0.5)^2) 

    if r <= 0.3
        return 9.0
    else 
        return 1.0
    end
end;
p = 1.0;
h = 1.0;
N = 100;

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
bz = collect(range(-π/p, π/p, N))

## Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=7)

## Plot the dispersion diagram
plot_diagram(bz, μ, period=p)