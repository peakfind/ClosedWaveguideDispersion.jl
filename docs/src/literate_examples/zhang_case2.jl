# # Inhomogeneous case

# ## Problem

# This problem is the second example in [zhang2021](@ref) page .
# + ``q(x_{1}, x_{2}) = 9`` in a disk with the center ``(0, 0.5)`` and radius ``0.3`` and ``q(x_{1}, x_{2}) = 1`` outside the disk.
# + Neumann boundary condition on ``\partial \Omega``.

# ## Code

# We need the same packages as the homogeneous case
using Ferrite: Lagrange, RefTriangle
using ClosedWaveguideDispersion

# Here we need to define a new refractive index. 
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

# Then we generate the mesh and set up variables for finite element method.
## Set up the grid 
grid = setup_grid(lc=0.05, period=p, height=h)

## Define the interpolation 
ip = Lagrange{RefTriangle, 1}()

## Set up the FE values: CellValues 
cv = setup_fevs(ip)

## Set up the DofHandler 
dh = setup_dofs(grid, ip)

# For the boundary condition, we only need to use [`setup_bdcs`](@ref) to impose the periodic boundary condition in ``x_{1}`` direction.
## Set the boundary conditions
cst = setup_bdcs(dh, period=p)

# We assemble the matrices for the generalized linear eigenvalue problems.
## Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

# Then we solve the generalized eigenvalue problems at each ``\alpha`` in ``bz``,
## Discretize the Brillouin zone
bz = collect(range(-π/p, π/p, N))

## Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=7)

# Finally, we can get our dispersion diagram.
## Plot the dispersion diagram
plot_diagram(bz, μ, period=p)