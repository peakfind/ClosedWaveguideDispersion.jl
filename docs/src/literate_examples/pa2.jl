# # PA-Two

# ## Problem

# PA-Two is an example in [ehrhardt-sun-zheng2009](@cite) and satisfies the following conditions:
# + a rectangular hole of size ``0.5 \times 0.5`` in the cell center;
# + Neumann boundary condition on ``\partial \Omega``;
# + homogeneous Dirichlet boundary condition at the hole boundary.


# ## Code
# For this example, we also need to load Gmsh and FerriteGmsh to generate the mesh
using Gmsh
using Ferrite
using FerriteGmsh
using ClosedWaveguideDispersion

# We need to customize our mesh and boundary conditions.
function my_grid(;lc=0.05, period=1.0, height=1.0, holewidth=0.5, holeheight=0.5)
    ## Initialize gmsh 
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    
    ## Add the points 
    p1 = gmsh.model.geo.addPoint(-period/2, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(period/2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(period/2, height, 0, lc)
    p4 = gmsh.model.geo.addPoint(-period/2, height, 0, lc)
    p5 = gmsh.model.geo.addPoint(-holewidth/2, (height - holeheight)/2, 0, lc)
    p6 = gmsh.model.geo.addPoint(holewidth/2, (height - holeheight)/2, 0, lc)
    p7 = gmsh.model.geo.addPoint(holewidth/2, (height + holeheight)/2, 0, lc)
    p8 = gmsh.model.geo.addPoint(-holewidth/2, (height + holeheight)/2, 0, lc)
    
    ## Add the lines 
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)
    
    ## Create loops and the domain 
    outerLoop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    innerLoop = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    domain = gmsh.model.geo.addPlaneSurface([outerLoop, innerLoop])
    
    ## Synchronize the model
    gmsh.model.geo.synchronize()

    ## Create the physical domains 
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(1, [l5, l6, l7, l8], -1, "Γ")
    gmsh.model.addPhysicalGroup(2, [domain], -1, "Ω")

    ## Set Periodic boundary condition
    transform = [1, 0, 0, period, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], transform) 

    ## Generate a 2D mesh
    gmsh.model.mesh.generate(2)

    grid = mktempdir() do dir 
        path = joinpath(dir, "hole.msh")
        gmsh.write(path)
        togrid(path)
    end
   
    ## Finalize the Gmsh library
    gmsh.finalize()
    
    return grid
end;

# We still use homogeneous medium in the waveguide.
function n(x)
    return 1.0
end;

# We need to implement a new function to impose the boundary condition.
function my_bdcs(dh::DofHandler; period=1.0)
    cst = ConstraintHandler(dh)
    
    ## Periodic boundary condition on the left and right 
    ## side of the periodic cell
    pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Ferrite.Vec{2}((period, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)
    
    ## Set Dirichlet boundary condition on the inner boundary
    dfacets = getfacetset(dh.grid, "Γ")
    dbc = Dirichlet(:u, dfacets, x -> 0)
    add!(cst, dbc)
    
    close!(cst)
end

# Set parameters
p = 1.0;
h = 1.0;
hw = 0.5;
hh = 0.5;
N = 150;

# Set up the grid 
grid = my_grid(lc=0.05, period=p, height=h, holewidth=hw, holeheight=hh)

# basic steps need by Ferrite.jl
## Define the interpolation 
ip = Lagrange{RefTriangle, 1}()
## Set up the FE values: CellValues 
cv = setup_fevs(ip)
## Set up the DofHandler 
dh = setup_dofs(grid, ip)

## Set the boundary conditions
cst = my_bdcs(dh, period=p)

## Allocate the matrices
A = allocate_matries(dh, cst)
B = allocate_matries(dh, cst)

## Discretize the Brillouin zone
bz = collect(range(-π/p, 2π/p, N))

## Calculate the dispersion diagram
μ = calc_diagram(cv, dh, cst, A, B, n, bz, nevs=5)

## Plot the dispersion diagram
plot_diagram(bz, μ, period=p)