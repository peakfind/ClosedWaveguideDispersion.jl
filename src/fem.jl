# Finite element method codes based on Ferrite.jl

"""
    setup_grid(;lc=0.05, period=1.0, height=1.0)

Generate the mesh for the periodic cell by using Gmsh.

# Arguments

- `lc`: the mesh size
- `period`: the period of the periodic closed waveguide
- `height`: the height of the periodic closed waveguide
"""
function setup_grid(;lc=0.05, period=1.0, height=1.0)
    # Initialize Gmsh 
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    
    # Add the points 
    p1 = gmsh.model.geo.addPoint(-period/2, 0, 0, lc)
    p2 = gmsh.model.geo.addPoint(period/2, 0, 0, lc)
    p3 = gmsh.model.geo.addPoint(period/2, height, 0, lc)
    p4 = gmsh.model.geo.addPoint(-period/2, height, 0, lc)

    # Add the lines 
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    
    # Create the loop and the surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Synchronize the model
    gmsh.model.geo.synchronize()
    
    # Create the physical groups 
    gmsh.model.addPhysicalGroup(1, [l1], -1, "bottom")
    gmsh.model.addPhysicalGroup(1, [l2], -1, "right")
    gmsh.model.addPhysicalGroup(1, [l3], -1, "top")
    gmsh.model.addPhysicalGroup(1, [l4], -1, "left")
    gmsh.model.addPhysicalGroup(2, [surf], -1, "Cell")

    # Set the periodic boundaries 
    transform = [1, 0, 0, period, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [l2], [l4], transform)

    # Generate a 2D mesh 
    gmsh.model.mesh.generate(2)
    
    # Read the .msh file by FerriteGmsh
    grid = mktempdir() do dir 
        path = joinpath(dir, "mesh.msh")
        gmsh.write(path)
        togrid(path)
    end

    # Finalize the Gmsh library 
    gmsh.finalize() 

    return grid
end

function setup_fevs(ip)
    qr = QuadratureRule{RefTriangle}(2)
    cv = CellValues(qr, ip)

    return cv
end

function setup_dofs(grid::Grid, ip)
    dh = DofHandler(grid)
    add!(dh, :u, ip)
    close!(dh)
    
    return dh
end

"""
    setup_bdcs(dh::DofHandler; period=1.0)

Impose the periodic boundary condition on the two vertical boundaries on the periodic cell.
"""
function setup_bdcs(dh::DofHandler; period=1.0)
    cst = ConstraintHandler(dh)
    
    # Periodic boundary condition on the left and right 
    # side of the periodic cell
    pfacets = collect_periodic_facets(dh.grid, "right", "left", x -> x + Ferrite.Vec{2}((period, 0.0)))
    pbc = PeriodicDirichlet(:u, pfacets)
    add!(cst, pbc)
    close!(cst)
    
    return cst
end

function allocate_matries(dh::DofHandler, cst::ConstraintHandler)
    sp = init_sparsity_pattern(dh)
    add_cell_entries!(sp, dh)
    add_constraint_entries!(sp, cst)
    K = allocate_matrix(SparseMatrixCSC{ComplexF64, Int}, sp)
    
    return K
end

"""
    assemble_A(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, α::Float64)

TBW
"""
function assemble_A(cv::CellValues, dh::DofHandler, A::SparseMatrixCSC, α::Float64)
    # Preallocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Ae = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler
    assembler = start_assemble(A)
    
    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell 
        reinit!(cv, cell)
        
        # Reset local matrix to 0.0 + 0.0im
        fill!(Ae, 0.0 + 0.0im)
        
        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)
                ∇v = shape_gradient(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    ∇u = shape_gradient(cv, qp, j)
                    
                    # Compute the local matrix according to the variational formulation
                    Ae[i, j] += (∇u ⋅ ∇v - (2im * α * ∇u[1] * v) + (α^2) * u * v) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Ae)
    end
    
    return A
end

"""
    assemble_B(cv::CellValues, dh::DofHandler, B::SparseMatrixCSC, n::Function)

TBW
"""
function assemble_B(cv::CellValues, dh::DofHandler, B::SparseMatrixCSC, n::Function)
    # Preallocate the local matrix
    n_basefuncs = getnbasefunctions(cv)
    Be = zeros(ComplexF64, n_basefuncs, n_basefuncs)
    
    # Create an assembler
    assembler = start_assemble(B)
    
    # Loop over all cells 
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell 
        reinit!(cv, cell)
        
        # Reset local matrix to 0.0 + 0.0im
        fill!(Be, 0.0 + 0.0im)
        
        # Get the coordinates of this cell 
        coords = getcoordinates(cell)

        # Loop over quadrature points
        for qp in 1:getnquadpoints(cv)
            dx = getdetJdV(cv, qp)
            
            # Get the coordinates of the quadrature point
            # and evaluate the refractive index at this point
            coords_qp = spatial_coordinate(cv, qp, coords)
            ri = n(coords_qp)
            
            # Loop over test shape functions 
            for i in 1:n_basefuncs
                v = shape_value(cv, qp, i)

                # Loop over trial shape functions 
                for j in 1:n_basefuncs
                    u = shape_value(cv, qp, j)
                    
                    # Compute the local matrix according to the variational formulation
                    Be[i, j] += (ri * u * v) * dx
                end
            end
        end
        
        assemble!(assembler, celldofs(cell), Be)
    end
    
    return B
end