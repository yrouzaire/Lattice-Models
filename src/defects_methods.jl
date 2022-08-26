include("../src/lattice_general.jl");
include("../src/models.jl");
include("../src/core_methods.jl");

function arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    Note that the inputs thetas are within [0,2π] =#
    dtheta = theta2 - theta1
    dtheta_abs = abs(theta2 - theta1)

    shortest_unsigned_arclength = min(symm-dtheta_abs,dtheta_abs)
    if dtheta_abs ≤ symm/2
        signe = sign(dtheta)
    else
        signe = -sign(dtheta)
    end
    return signe*shortest_unsigned_arclength
end

function get_vorticity(thetasmodpi::Matrix{T},i::Int,j::Int,L::Int,sym)::T where T<:AbstractFloat
    #= Note : thetasmodpi = mod.(thetas,symm)
        By convention, the neighbours have the following ordering
           3   2
        4    x   1
           5   6
        =#

    angles_corners = get_neighbours_triangular(thetasmodpi,L,i,j,is_in_bulk(i,j,L))
    perimeter_covered = 0.0
    for i in 1:length(angles_corners)-1
        perimeter_covered += arclength(angles_corners[i],angles_corners[i+1],sym)
    end
    perimeter_covered += arclength(angles_corners[end],angles_corners[1],sym)
    charge = round(perimeter_covered/2π,digits=1)
    return charge
end

function spot_defects(thetas::Matrix{Float32},T::Number,BC::String)
function spot_defects(model::AbstractModel,lattice::AbstractLattice)
    L = model.L
    lattice.periodic ? range_bc = 1:L : range_bc = 2:L-1
    list_vortices_plus  = Tuple{Int,Int}[]
    list_vortices_minus = Tuple{Int,Int}[]

    #= Relaxation of a fictitious system (it won't evolve the actual system)
    at 0 temperature so that vortices are easier to detect (no spurious doubles because of temperature variations)
    After some crude benchamarks, it seems that the number of relaxation loops requiered to erase all the doubles
    is approximately 20T. Though, if possible, the lesser the better because vortices might move during this process -> inadéquation
    entre la position affichée des défauts et le heatmap.  =#
    # thetas_prerelaxed = remove_isolate(copy(thetas),L)
    # for n in 1:10 thetas_prerelaxed = align(thetas_prerelaxed,L,0,0,"polar",BC) end
    # thetasmodpi = mod.(fill_relax_holes(thetas_prerelaxed),π)
    # thetasmodpi = cg_gaussian_fiecld(remove_isolate(copy(thetas),L))
    thetasmodpi = mod.(model.thetas,sym(model))
    for i in range_bc
        for j in range_bc
            q = get_vorticity(thetasmodpi,i,j,L)
            if     q > + 0.1 push!(list_vortices_plus,(i,j)) # we want to keep ±½ defects, and not rounding errors
            elseif q < - 0.1 push!(list_vortices_minus,(i,j))
            end
        end
    end
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#

    elements_to_keep_plus  = trues(length(list_vortices_plus))
    elements_to_keep_minus = trues(length(list_vortices_minus))
    seuil = 4
    for i in 2:length(list_vortices_minus)
        for j in 1:i-1
            if dist(list_vortices_minus[i],list_vortices_minus[j],L) ≤ seuil
                elements_to_keep_minus[i] = 0
                # break # will only break the inner for loop and continue the outer one
            end
        end
    end
    for i in 2:length(list_vortices_plus)
        for j in 1:i-1
            if dist(list_vortices_plus[i],list_vortices_plus[j],L) ≤ seuil
                elements_to_keep_plus[i] = 0
                # break # will only break the inner for loop and continue the outer one
            end
        end
    end
    # return list_vortices_plus,list_vortices_minus
    return list_vortices_plus[elements_to_keep_plus],list_vortices_minus[elements_to_keep_minus]
end

function number_defects(thetas,T,BC)
    a,b = spot_defects(thetas,T,BC)
    return length(a) + length(b)
end

function theta_mid(x::T,y::T,symm)::T  where T<:AbstractFloat
    arcl = arclength(x,y,symm)
    return x + arcl/2 , y - arcl/2
end
