"""
    SquareLattice{T}

2D Square lattice described by x and y.
```   
(0, y)
     ^ 
     |
     |
     |
     ---------> (x, 0) 
```
"""
struct SquareLattice{T}
    x::T
    y::T
end

"""
    IrreducibleBrillouin{T}

Irreducible Brillouin zone of 2D square lattices. It is a triangle with three 
vertices: Γ, X, and M.

# Fields

- `Γ`: Center of the Brillouin zone ``(0, 0)``
- `X`: ``(\\pi / x, 0)``
- `M`: ``(\\pi / x, \\pi / y)``
"""
struct IrreducibleBrillouin{T}
    Γ::SVector{2, T}
    X::SVector{2, T}
    M::SVector{2, T}
end

function IrreducibleBrillouin(sl::SquareLattice{T}) where {T}
    π_x = T(π) / sl.x
    π_y = T(π) / sl.y

    Γ = @SVector zeros(T, 2)
    X = @SVector [π_x, zero(T)]
    M = @SVector [π_x, π_y]

    return IrreducibleBrillouin(Γ, X, M)
end

"""
    get_discrete_irreducibleBrillouin(ibz::IrreducibleBrillouin{T}, n::Int64) where{T}

Get discrete points on the boundary of the irreducible Brillouin zone (a 
triangle for square lattices): Γ -> X -> M -> Γ.

                 M
                /|   
               / |
             /   |
       d   /     | π_y
         /       |
       /         |
     /           |
    --------------
    Γ     π_x    X

# Arguments

- `ibz`: the irreducible Brillouin zone
- `n`: the number of points in the interior of each edge

# Output

- `dibz`: the points along the boundary of the irreducible Brillouin zone
- `para`: a real number in ``[0, 1]``. which is the parameter corresponding 
          to the points in `dibz`. We need `para` when we plot the band structure
"""
function get_discrete_irreducibleBrillouin(ibz::IrreducibleBrillouin{T}, n::Int64) where {T}
    n > 0 || throw(ArgumentError("n must be positive"))

    # Total number of the points along the Γ -> X -> M -> Γ
    N = 3 * n + 4

    # Preallocate
    dibz = Vector{SVector{2, T}}(undef, N)
    para = Vector{T}(undef, N)

    # Precompute some constants
    π_x = ibz.X[1]
    π_y = ibz.M[2]
    hx = π_x / (n + 1)
    hy = π_y / (n + 1)
    d = StaticArrays.norm(ibz.M)

    # Perimeter of the Irreducible Brillouin zone
    l = π_x + π_y + d

    # Γ -> X (except X)
    dibz[1] = ibz.Γ
    para[1] = zero(T)

    for i in 2:(n + 1)
        tx = (i - 1) * hx
        p = SVector{2, T}(tx, zero(T))
        dibz[i] = p
        para[i] = tx / l
    end

    # X -> M (except M)
    dibz[n + 2] = ibz.X
    para[n + 2] = π_x / l

    for i in (n + 3):(2 * n + 2)
        ty = (i - n - 2) * hy
        p = SVector{2, T}(π_x, ty)
        dibz[i] = p
        para[i] = (ty + π_x) / l
    end

    # M -> Γ
    dibz[2 * n + 3] = ibz.M
    para[2 * n + 3] = (π_x + π_y) / l

    for i in (2 * n + 4):(3 * n + 3)
        tx = (N - i) * hx
        ty = (N - i) * hy
        td = (i - 2 * n - 3) * d / (n + 1)
        p = SVector{2, T}(tx, ty)
        dibz[i] = p
        para[i] = (π_x + π_y + td) / l
    end

    dibz[N] = ibz.Γ
    para[N] = one(T)

    return dibz, para
end
