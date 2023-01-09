## ------------------------ Lattices ------------------------
abstract type AbstractLattice end
abstract type Abstract1DLattice <: AbstractLattice end
abstract type Abstract2DLattice <: AbstractLattice end

mutable struct Chain1D <: Abstract1DLattice
    L::Int
    periodic::Bool
end
Chain1D(L::Int;periodic::Bool=true) = Chain1D(L,periodic)

mutable struct TriangularLattice <: Abstract2DLattice
    L::Int
    periodic::Bool
    single::Bool
    metric::String
end
TriangularLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = TriangularLattice(L,periodic,single,metric)

mutable struct SquareLattice <: Abstract2DLattice
    L::Int
    periodic::Bool
    single::Bool
    metric::String
end
SquareLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = SquareLattice(L,periodic,single,metric)
function number_nearest_neighbours(lattice::Abstract2DLattice)
    if isa(lattice,TriangularLattice) nnn = 6
    elseif isa(lattice,SquareLattice) nnn = 4
    end
    return nnn
end

## ------------------------ Functions ------------------------
function dist(lattice::Abstract2DLattice,pos1,pos2)
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


function distance_matrix(new,old,lattice::AbstractLattice)
    m_new,m_old = length(new),length(old)
    distance_matrix = zeros(m_new,m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i,j] = dist(lattice,new[i],old[j])
        end
    end
    return distance_matrix
end

#= Pour les offsets, on commence par le voisin du bas, dans le référentiel
de la matrice, ie. axe j vers la droite et axe i vers le bas
La raison: parce qu'apres un heatmap(thetas'), equivalent à
une rotation de 90° antihoraire, le premier element du vecteur "offset" sera
à droite, cohérent avec l'intuition qu'on s'en fait pour \theta = 0
(qui sera donc projeté sur le premier element de offsets) =#
function offsets(lattice::TriangularLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if even return [(1,0) , (1,1) , (0,1)  , (-1,1)  , (-1,0) ,  (0,-1) ]
    else    return [(1,0) , (0,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1) ]
    end
end

function offsets(lattice::SquareLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if lattice.metric in ["manhattan","euclidian"]
        return [(1,0), (0,1) , (-1,0) , (0,-1) ]
    elseif lattice.metric =="chebychev"
        return [(1,0) , (1,1) , (0,1) , (-1,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1)]
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

true_position(lattice::SquareLattice,i,j) = (i,j)
function true_position(lattice::TriangularLattice,i,j)
    if iseven(i) return (i,j+0.5)
    else return (i,j)
    end
end

function distance_to_border(thetas::Matrix{<:AbstractFloat},i,j)
    L = size(thetas,1)
    distance_to_left   = i-1
    distance_to_right  = L-i
    distance_to_bottom = j-1
    distance_to_top    = L-j
    return minimum([distance_to_top,distance_to_left,distance_to_right,distance_to_bottom])
end

function add_2_positions(pos1::Tuple{T,T},pos2::Tuple{T,T},L::T,should_take_mod::Bool)::Tuple{T,T} where T<:Int
    if should_take_mod return mod1.(pos1 .+ pos2 ,L)
    else return pos1 .+ pos2
    end
end

function mean_2_positions(pos1,pos2,L,should_take_mod::Bool=true)
    a,b = pos1 ; x,y = pos2

    dx = (x - a) #; dx = min(dx,L-dx)
    dy = (y - b) #; dy = min(dy,L-dy)

    if should_take_mod
        if abs(L-dx) < abs(dx) dx = -(L-dx) end
        # if L-dx < dx dx = -(L-dx) end
        if abs(L-dy) < abs(dy) dy = -(L-dy) end
        # if L-dy < dy dy = -(L-dy) end
        estimate_loc_collision = mod1.((a,b) .+ 0.5.*(dx,dy),L)
        # estimate_loc_collision = mod1.(0.5.*(a,b) .+ 0.5.*(x,y),L)
        return estimate_loc_collision
    end
    return (a,b) .+ 0.5.*(dx,dy) # if should_take_mod == false
end
# l = 100
# mean_2_positions((50,50),(60,60),l) == (55,55)
# mean_2_positions((10,10),(90,90),l) == (100,100)
# mean_2_positions((49,66),(51,61),l) == (50.0, 63.5)

function mean_N_positions(vec_pos,L,should_take_mod::Bool=true)
    averaged_pos = vec_pos[1]
    for i in 2:length(vec_pos)
        averaged_pos = mean_2_positions(averaged_pos,vec_pos[i],L,should_take_mod)
    end
    return averaged_pos
end

# coordinate-free methods based on scalar products (does not force SquareLattice)
function get_div_rot(thetas::Matrix{T},lattice::Abstract2DLattice) where T<:AbstractFloat # coordinate-free
    L = size(thetas,1)
    divergence = NaN*zeros(L,L)
    rotational = NaN*zeros(L,L)
    for j in 1:L , i in 1:L
        divergence[i,j] , rotational[i,j] = get_div_rot(thetas,i,j,lattice)
    end
    return divergence , rotational
end

function get_div_rot(thetas::Matrix{T},i,j,lattice::Abstract2DLattice) where T<:AbstractFloat # coordinate-free
    L = size(thetas,1)

    dummy_rho = one(T) # so that the NaN do not get filtered out, I need them here
    rho1_model = LangevinXY{T}(zero(T),"polar",zero(T),zero(T),dummy_rho)

    if isa(lattice,TriangularLattice) cst = π/3 # ; surface_unit_cell = 3sqrt(3)/2
    elseif isa(lattice,SquareLattice) cst = π/2 # ; surface_unit_cell = 2
    end

    angles_neighbours = get_neighbours(thetas,rho1_model,lattice,i,j,is_in_bulk(i,j,L))
    div = 0.0
    rot = 0.0
    for k in 1:length(angles_neighbours)
        if !isnan(angles_neighbours[k]) # normally it should not be, since preconditionning! first, but one never knows.
            div += cos((k-1)*cst) * cos(angles_neighbours[k]) + sin((k-1)*cst) * sin(angles_neighbours[k])
            rot += cos(k*cst) * cos(angles_neighbours[k]) + sin(k*cst) * sin(angles_neighbours[k])
        end
    end
    return div,rot
    #= No division by surface_unit_cell because eventually
    we do not care about the actual numerical value,
    only its sign. =#
end

# methods based on derivatives (forces SquareLattice)
# function get_div_rot(thetas::Matrix{T}) where T<:AbstractFloat # derivative based
#     L = size(thetas,1)
#     divergence = NaN*zeros(L,L)
#     rotational = NaN*zeros(L,L)
#     for j in 1:L , i in 1:L
#         divergence[i,j] , rotational[i,j] = get_div_rot(thetas,i,j)
#     end
#     return divergence , rotational
# end
#
# function get_div_rot(thetas::Matrix{T},i,j) where T<:AbstractFloat # derivative based
#     L = size(thetas,1)
#
#     dummy_rho = one(T) # so that the NaN do not get filtered out, I need them here
#     rho1_model = XY{T}(zero(T),"polar",zero(T),zero(T),dummy_rho)
#     square_lattice = SquareLattice(L,periodic=false)
#
#     angles_neighbours = get_neighbours(thetas,rho1_model,square_lattice,i,j,is_in_bulk(i,j,L))
#     # TODO check the order of angles_neighbours[1 or 3] and whether 2pi systematiquement
#     dx_theta = arclength(angles_neighbours[3],angles_neighbours[1],2pi)/2
#     dy_theta = arclength(angles_neighbours[4],angles_neighbours[2],2pi)/2
#     theta = thetas[i,j]
#     div = -sin(theta)*dx_theta + cos(theta)*dy_theta
#     rot =  cos(theta)*dx_theta + sin(theta)*dy_theta
#     return div,rot
# end

function zoom(thetas::Matrix{<:Number},lattice::Abstract2DLattice,i,j,window=WINDOW)
    #= Returns a 2window+1 portion of thetas (centered on i,j),
    together with a bool to say whether the operation has worked =#
    L = lattice.L
    i = round(Int,i)
    j = round(Int,j)

    if at_least_at_distance_X_from_border(i,j,L,X=window+1)
        thetas_zoom = thetas[i-window:i+window,j-window:j+window]
        no_problem_go_ahead = true
    else
        if lattice.periodic
            shift = 20
            thetas_shifted = copy(thetas)
            for i in 1:L thetas_shifted[i,:] = circshift(thetas_shifted[i,:],shift) end
            for i in 1:L thetas_shifted[:,i] = circshift(thetas_shifted[:,i],shift) end
            return zoom(thetas_shifted,lattice,mod1(i+shift,L),mod1(j+shift,L),window)

        else
            no_problem_go_ahead = false
            thetas_zoom = NaN
        end
    end
    return no_problem_go_ahead,thetas_zoom
end

## Rotations
rotate_clockwise90(thetas::Matrix{<:Number}) = rotr90(thetas .- pi/2)
rotate_counterclockwise90(thetas::Matrix{<:Number}) = rotl90(thetas .+ pi/2)
rotate_180(thetas::Matrix{<:Number}) = rot180(thetas .+ pi)
function randomly_rotate(thetas::Matrix{T})::Matrix{T} where T<:Number
    u = rand()
    if     u < 0.25  return thetas
    elseif u < 0.50  return rotate_clockwise90(thetas)
    elseif u < 0.75  return rotate_counterclockwise90(thetas)
    else             return rotate_180(thetas)
    end
end
