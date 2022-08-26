export TriangularLattice

abstract type AbstractLattice end
struct TriangularLattice <: AbstractLattice
    L::Int
    periodic::Bool
    single::Bool
    metric::String
end
TriangularLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = TriangularLattice(L,periodic,single,metric)

function dist(lattice::AbstractLattice,pos1::Tuple{T,T},pos2::Tuple{T,T}) where T<:Number
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    if lattice.periodic
        dx = min(dx,lattice.L-dx)
        dy = min(dy,lattice.L-dy)
    end
    if lattice.metric == "euclidian"     return sqrt(dx^2 + dy^2)
    elseif lattice.metric == "manhattan" return dx + dy
    end
end

function linear_to_square_index(n::Int,L::Int)
    i = div(n-1,L) + 1  # 1 ≤ i ≤ L
    j = mod1(n,L)       # 1 ≤ j ≤ L
    #= formula for reversing the linear indexing. i is the quotient and j
    the reminder of the euclidian division of n by L and the corrections
    deal with the 1-indexing of Julia =#
    return i,j
end

function square_to_linear_index(i::Int,j::Int,L::Int)::Int
    return L*(i-1) + j # 1 ≤ n ≤ L^2
end

function closetoborders(i::Int,j::Int,L::Int,c::Int=5)
    if     i ≤ c     return true
    elseif j ≤ c     return true
    elseif i + c ≥ L return true
    elseif j + c ≥ L return true
    else return false
    end
end

function ontheborder(i::Int,j::Int,L::Int)::Bool
    if i == 1 || j == 1 || i == L || j == L
        return true
    else
        return false
    end
end

is_in_bulk(i::Int,j::Int,L::Int) = !ontheborder(i::Int,j::Int,L::Int)
