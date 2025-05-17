"""
    calc_diagram(cv::CellValues, dh::DofHandler, cst::ConstraintHandler, A::SparseMatrixCSC, B::SparseMatrixCSC, n::Function, bz; nevs::Int=6)

TBW
"""
function calc_diagram(cv::CellValues, dh::DofHandler, cst::ConstraintHandler, A::SparseMatrixCSC, B::SparseMatrixCSC, n::Function, bz; nevs::Int=6)
    m = length(bz)
    μ = zeros(m, nevs)
    
    # Assemble the matrix B, since B is independent of α
    B = assemble_B(cv, dh, B, n)
    
    # Impose the boundary conditions
    apply!(B, cst)
    
    # Assemble the matrix A at different α in bz
    for (i, α) in enumerate(bz) 
        # Assemble A at α 
        A = assemble_A(cv, dh, A, α)

        # Impose the boundary conditions
        apply!(A, cst)
        
        # Solve the generalized eigenvalue problem by Arpack.jl
        λ, _ = eigs(A, B, nev=nevs, which=:SR) 
        λ = real(λ)
        μ[i, :] = λ
        
        # Reset A to 0.0 + 0.0im
        fill!(A, 0.0 + 0.0im)
    end
    
    return μ
end