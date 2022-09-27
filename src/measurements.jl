include("lattices.jl");
include("models.jl");

function OP(thetas::Matrix{T}) where T<:AbstractFloat
    tmp_polar   = Complex(0)
    tmp_nematic = Complex(0)
    for theta in thetas
        if !isnan(theta)
            tmp_polar   += exp(im*theta)
            tmp_nematic += exp(2im*theta)
        end
    end
    return abs(tmp_polar/length(thetas)) , abs(tmp_nematic/length(thetas))
end

function local_mag(thetas::Matrix{T},L::Int,l::Int)::Matrix{T} where T<:AbstractFloat
    # L is the size of the system
    # l is the size of the local box for averaging the magnetisation
    @assert isinteger(L/l)
    c = Int(L/l)
    δ = zeros(2,c^2)
    for j in 1:c
        for i in 1:c
            δ[:,square_to_linear_index(i,j,c)] = get_Delta(thetas[1+l*(i-1):i*l, 1+l*(j-1):j*l])
        end
    end
    return δ
end

function corr(thetas::Matrix{T},model::AbstractModel,lattice::Abstract2DLattice) where T<:AbstractFloat
    L = lattice.L ; Lover2 = round(Int,L/2,RoundDown)
    matrix_corr = Matrix{T}(undef,Lover2,L^2)

    if lattice.periodic nmax = Lover2
    else nmax = distance_to_border(thetas,i,j)-1
    end

    for m in 1:L^2
        i,j = linear_to_square_index(m,L)
        for n in 1:nmax
            matrix_corr[n,m] = corr(thetas,model,lattice,i,j,n)
        end
    end
    return nanmean(matrix_corr,2)[:,1]
end

function corr(thetas::Matrix{T},model::AbstractModel,lattice::Abstract2DLattice,i::Int,j::Int,n::Int) where T<:AbstractFloat
    angle_neighbours_at_distance_n = [thetas[mod1(i+n,L),j],
                                      thetas[mod1(i-n,L),j],
                                      thetas[i,mod1(j+n,L)],
                                      thetas[i,mod1(j-n,L)] ]
    filter!(!isnan,angle_neighbours_at_distance_n)
    model.symmetry == "polar" ? symm = 1 : symm = 2
    if length(angle_neighbours_at_distance_n) > 0
        return mean(cos,symm*(angle_neighbours_at_distance_n .- thetas[i,j]))
    else
        return NaN
    end
end

# Method with correct_sites way too slow, so let's go back to a flawed but tractable solution : use the same algo for TriangularLattice and SquareLattice
# function corr(thetas::Matrix{T},model::AbstractModel,lattice::TriangularLattice,i::Int,j::Int,n::Int) where T<:AbstractFloat
#     # correct_sites = []
#     # for a in -n:+n , b in -n:+n # scann the square with chebychev distance ≤ n
#     #     first_selection = abs(a) + abs(b) ≥ n # to exclude already some inner part of the square
#     #     correct_distance = (dist(lattice,true_position(lattice,a,b),true_position(lattice,i,j)) == n)
#     #     if first_selection && correct_distance
#     #         push!(correct_sites,(a,b))
#     #     end
#     # end
#     angle_neighbours_at_distance_n = [thetas[mod1(i+n,L),j],
#                                       thetas[mod1(i-n,L),j],
#                                       thetas[i,mod1(j+n,L)],
#                                       thetas[i,mod1(j-n,L)],
#                                       thetas[mod1(i-n,L),mod1(j-n,L)],
#                                       thetas[mod1(i+n,L),mod1(j+n,L)],
#                                       thetas[mod1(i-n,L),mod1(j+n,L)],
#                                       thetas[mod1(i+n,L),mod1(j-n,L)] ]
#     posij = true_position(lattice,i,j)
#     filter!(x->(dist(lattice,true_position(lattice,x...),posij) == n),angle_neighbours_at_distance_n)
#     filter!(!isnan,angle_neighbours_at_distance_n)
#     model.symmetry == "polar" ? symm = 1 : symm = 2
#     if length(angle_neighbours_at_distance_n) > 0
#         return mean(cos,symm*(angle_neighbours_at_distance_n .- thetas[i,j]))
#     else
#         return NaN
#     end
# end

function corr_length(C::Vector{T},rs=1:length(C);seuil=exp(-1))::T where T<:AbstractFloat # from a time series, returns the correlation length ()
    i_after = findfirst(x->x<seuil,C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
    else
    ξ = NaN
    end
    return ξ
end
