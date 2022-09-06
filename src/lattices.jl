## ------------------------ Lattices ------------------------
abstract type AbstractLattice end

mutable struct TriangularLattice <: AbstractLattice
    const L::Int
    periodic::Bool
    const single::Bool
    const metric::String
end
TriangularLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = TriangularLattice(L,periodic,single,metric)

mutable struct SquareLattice <: AbstractLattice
    const L::Int
    periodic::Bool
    const single::Bool
    const metric::String
end
SquareLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = SquareLattice(L,periodic,single,metric)

## ------------------------ Functions ------------------------
function dist(lattice::AbstractLattice,pos1::Tuple{T,T},pos2::Tuple{T,T}) where T<:Number
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    if lattice.periodic
        dx = min(dx,lattice.L-dx)
        dy = min(dy,lattice.L-dy)
    end
    if     lattice.metric == "euclidian" return sqrt(dx^2 + dy^2)
    elseif lattice.metric == "manhattan" return dx + dy
    elseif lattice.metric == "chebychev" return max(dx,dy)
    else error("Unknown metric !")
    end
end

function offsets(lattice::TriangularLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if even return [(0,1) , (-1,1) , (-1,0) , (0,-1)  , (1,0)  , (1,1)]
    else    return [(0,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1) , (1,0)]
    end
end

function offsets(lattice::SquareLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if lattice.metric in ["manhattan","euclidian"]
        return [(0,1) , (-1,0) , (0,-1) , (1,0)]
    elseif lattice.metric =="chebychev"
        return [(0,1) , (-1,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1), (1,0) , (1,1)]
    else error("Unknown metric !")
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

function ontheborder(i::Int,j::Int,L::Int)::Bool
    if i == 1 || j == 1 || i == L || j == L
        return true
    else
        return false
    end
end

is_in_bulk(i::Int,j::Int,L::Int)::Bool = !ontheborder(i,j,L)
function at_least_at_distance_X_from_border(i::Int,j::Int,L::Int;X=2)::Bool
    if i < X || j < X || i > L-X+1 || j > L-X+1
        return false
    else
        return true
    end
end

function add_2_positions(pos1::Tuple{T,T},pos2::Tuple{T,T},L::T,should_take_mod::Bool)::Tuple{T,T} where T<:Int
    if should_take_mod return mod1.(pos1 .+ pos2 ,L)
    else return pos1 .+ pos2
    end
end

# function add_2_positions(pos1::CartesianIndex{2},pos2::CartesianIndex{2},L::Int,should_take_mod::Bool)::CartesianIndex{2}
#     if should_take_mod return mod1.((pos1 + pos2).I ,L)
#     else return pos1 + pos2
#     end
# end
#
# a = CartesianIndex(1,1)
# b = CartesianIndex(1,10)
# mod.((a+b).I,10)
