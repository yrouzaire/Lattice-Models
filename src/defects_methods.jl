include("lattices.jl");
include("models.jl");
include("core_methods.jl");

using Hungarian
import Random.shuffle

function arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING
    Note that the inputs thetas need to lie within [0,π] or [0,2π], depending on the symmetry of the model =#
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

function get_vorticity(thetasmod::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice,i::Int,j::Int)::T where T<:AbstractFloat
    #= Note : thetasmod = mod.(thetas,symm*pi)
        By convention, for a SquareLattice, the neighbours have the following ordering
           2
        3  x  1
           4
       By convention, for a TriangularLattice, the neighbours have the following ordering
          3  2
        4  x  1
          5  6
        =#
    if     model.symmetry == "nematic" symm = T(pi)
    elseif model.symmetry == "polar"   symm = T(2pi)
    end
    angles_corners = invoke(get_neighbours,Tuple{Matrix{T},AbstractModel{T},typeof(lattice),Int,Int,Bool},   thetasmod,model,lattice,i,j,is_in_bulk(i,j,lattice.L))
    # the invoke trick above is so that the neighbors are all considered, even in the Non Reciprocal case
    perimeter_covered = 0.0
    for i in 1:length(angles_corners)-1
        perimeter_covered += arclength(angles_corners[i],angles_corners[i+1],symm)
    end
    if !isempty(angles_corners)
        perimeter_covered += arclength(angles_corners[end],angles_corners[1],symm)
        charge = round(perimeter_covered/2π,digits=1) # here 2π no matter the symmetry of the model. if nematic, will return half-integer defects
    else charge = NaN
    end
    return charge
end

function spot_defects(thetas::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice; seuil=4) where T<:AbstractFloat
    L = lattice.L
    if lattice.periodic range_bc = 1:L else range_bc = 2:L-1 end
    vortices_plus  = Tuple{Int,Int,T}[]
    vortices_minus = Tuple{Int,Int,T}[]

    if     model.symmetry == "nematic" modd = T(pi)
    elseif model.symmetry == "polar"   modd = T(2pi)
    end
    thetasmod = mod.(copy(thetas),modd)
    if model.rho < 1 preconditionning!(thetasmod,model,lattice) end
    for i in range_bc
        for j in range_bc
            q = get_vorticity(thetasmod,model,lattice,i,j)
            if     q > + 0.1 push!(vortices_plus,(i,j,q)) # we want to keep ±½ defects, and not rounding errors
            elseif q < - 0.1 push!(vortices_minus,(i,j,q))
            end
        end
    end
    return merge_duplicates(vortices_plus,lattice),merge_duplicates(vortices_minus,lattice)

    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    # vortices_to_keep_plus  = vortices_to_keep(vortices_plus,lattice,seuil)
    # vortices_to_keep_minus = vortices_to_keep(vortices_minus,lattice,seuil)
    #
    # return vortices_plus,vortices_minus
    # return vortices_plus[vortices_to_keep_plus],vortices_minus[vortices_to_keep_minus]
end

function merge_duplicates(list,lattice;radius=2)
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    pos    = [list[i][1:2] for i in each(list)]
    charge = [list[i][3]   for i in each(list)]
    dealt_with = falses(length(pos))

    merged_duplicates = []
    for i in 1:length(pos)
        if !dealt_with[i]
            tmp = []
            for j in i:length(pos) # includes defect "i"
                if dist(lattice,pos[i],pos[j]) ≤ radius && charge[i] == charge[j]
                    dealt_with[j] = true
                    push!(tmp,pos[j])
                end
            end
            push!(merged_duplicates,mean_N_positions(tmp,lattice.L,true))
        end
    end
    return merged_duplicates
end

function preconditionning!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::AbstractLattice)
    remove_isolates!(thetas,model,lattice)
    fill_holes!(thetas,model,lattice)
    relax!(thetas,model)
end

function remove_isolates!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::AbstractLattice)
    L = lattice.L
    for j in 1:L
        for i in 1:L
            if isempty(get_neighbours(thetas,model,lattice,i,j,is_in_bulk(i,j,L)))
                thetas[i,j] = NaN
            end
        end
    end
    return thetas
end

function fill_holes!(thetas::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice) where T<:AbstractFloat
    L = lattice.L
    # model.symmetry == "polar" ? symm = 2 : symm = 1
    while count(isnan,thetas) > 0
        for j in shuffle(1:L) , i in shuffle(1:L)
            if isnan(thetas[i,j])
                angle_neighbours = invoke(get_neighbours,Tuple{Matrix{T},AbstractModel{T},typeof(lattice),Int,Int,Bool},  thetas,model,lattice,i,j,is_in_bulk(i,j,L))
                # thetas[i,j] = angle(sum(exp,symm*im*angle_neighbours, init=0))/symm # moyenne des voisins
                if !isempty(angle_neighbours) thetas[i,j] = rand(angle_neighbours) end
            end
        end
    end
    return thetas
end



function vortices_to_keep(vortices,lattice,seuil)
    elements_to_keep = trues(length(vortices))
    for i in 2:length(vortices)
        for j in 1:i-1
            if dist(lattice,vortices[i][1:2],vortices[j][1:2]) ≤ seuil
                elements_to_keep[i] = false
                # break # will only break the inner for loop and continue the outer one
            end
        end
    end
    return elements_to_keep
end

# function spot_single_default_global(thetas::Matrix{T},model::AbstractModel,lattice::AbstractLattice)::Tuple{Tuple{Int16,Int16},T} where T<:AbstractFloat
#     L = size(thetas)[1]
#     list_defaults = Tuple{Int16,Int16,T}[]
#     for i in 2:L-1
#         for j in 2:L-1
#             q = get_vorticity(thetas,model,lattice,i,j)
#             if abs(q) > 0.1 # to avoid rounding errors
#                 push!(list_defaults,(i,j,q))
#             end
#         end
#     end
#     @assert length(list_defaults) == 1
#     return list_defaults[1][1:2],list_defaults[1][3]
# end


#is it even useful now ?
# function spot_single_default_local(thetas::Matrix{<:AbstractFloat},last_loc::T,known_loc::T,Q::Number,margin::Int=6)::Tuple{T,Bool} where T<:Tuple{Int16,Int16}
#     # V    is the location (x,y,q) of the default
#     #= The search for the new location of the default will scan a square of
#     dimension 2margin x 2margin (within the limits of the LxL lattice) ,
#     with the last location as center. The whole method is O(margin²). =#
#     L = size(thetas)[1]
#     positions = []
#
#
#     # Dimensions of the square
#     j_low,j_high = max(1,last_loc[2]-margin) , min(L,last_loc[2]+margin)
#     i_low,i_high = max(1,last_loc[1]-margin) , min(L,last_loc[1]+margin)
#     for j in j_low:j_high
#         for i in i_low:i_high
#             q = get_vorticity(thetas,i,j,L)
#             if abs(q) > 0.1
#                 push!(positions,(i,j,q))
#             end
#         end
#     end
#
#
#     #= In this list, there might be doubles/triples (2/3 locations for the
#     same physical vortex). We thus seek for numerically identified vortices
#     which are neighbours and with the same charge to delete them. =#
#     elements_to_keep = BitVector(ones(length(positions)))
#     for i in 2:length(positions)
#         for j in 1:i-1
#             if positions[i][3] == positions[j][3] && isneighbour(positions[i][1:2],positions[j][1:2],L)
#                 # if same charge and isneighbour(vortex_i,vortex_j) == true (for at least one j), we have a double, so we don't keep it
#                 elements_to_keep[i] = 0
#                 break # will only break the inner for loop and continue the outer one
#             end
#         end
#     end
#     positions = positions[elements_to_keep]
#     ℓ = length(positions)
#     # println("ℓ = $positions")
#
#     #= Summary of what follows :
#     ℓ = 0 : (i)   Controled Annihilation (encounter of the vortex with its antivortex)
#             (ii)  Unexpected Annihilation (encounter of the vortex with another antivortex
#             (iii) The (single) vortex left the lattice.
#             In all three case, throw an explicit error and treat it later on.
#
#     ℓ = 1 : All good ! One is pretty sure that the only detected vortex is the same
#     than the one of the previous timestep. The only case where we could be mistaken
#     is if a pair of vortex enters the square and that our previous vortex annihilates
#     with the newcomer antivortex. In this case, the only default remaining would be
#     the newcomer vortex, and we would be clueless about it. The only signature would be a
#     possible jump in space in the r(t) graph.
#
#     ℓ = 2 : If the other default is known/authorized, meaning that it is the former
#     antivortex of the default we currently work on, that's not important and we keep
#     going as if the vortex was alone. If the other defaut is NOT known/authorized,
#     it will perturb the displacment of our vortex : don't save the current position
#     by signalling it (alone = false)
#
#     ℓ ≥ 3 : We are sure to be in the bad case #2.  =#
#     if ℓ == 1
#         alone = true # no other default in the searched region
#         most_probable = positions[1][1:2]
#         #= Actually,  =#
#     elseif ℓ > 1
#         alone = false # by default. We deal we the special case ℓ = 2 just below.
#
#         if ℓ == 2 && ((known_loc[1],known_loc[2],-Q) in positions) ; alone = true ; end
#         #= If ℓ == 2, one has 2 possibilities :
#             - either the extra default is "known" (its the antidefault of
#             the one we currently consider), in which case turn alone = true.
#             - or it's not, in which case we leave alone = false.
#             In any case, the following routine to determine the location of
#             the default currenlty tracked is still completely valid.
#         Note that when using this function for the study of a single vortex,
#         one needs to provide an impossible known_loc, such as (-1,-1). =#
#
#         distances_to_last = [dist(positions[i][1:2],last_loc,L) for i in 1:ℓ]
#         # Kill candidates with opposite polarity
#         for i in 1:ℓ
#             element = positions[i]
#             if element[3] ≠ Q distances_to_last[i] = Inf end
#         end
#         most_probable = positions[sortperm(distances_to_last)][1]
#         #= Returns the position with the smallest distance to the last location,
#         hence we choose that one as the most probable candidate for the vortex we
#         are considering. =#
#
#     else # dealing with the case ℓ = 0
#
#         close_from_boundary = last_loc[1] < 4 || last_loc[1] > L - 4 || last_loc[2] < 4 || last_loc[2] > L - 4
#         if close_from_boundary                      error("The vortex left the lattice.")
#         elseif dist(last_loc,known_loc,L) ≤ sqrt(5) error("Controlled annihilation.") # sqrt(5) is arbitrary
#         else                                        error("Unexpected annihilation.") end
#     end
#     # println(1)
#     # println(typeof(Int16.(most_probable[1:2])))
#     return most_probable[1:2],alone
# end

## Defect Tracking (with a lot of defects, creation and annihilation)
mutable struct Defect
    id::Int
    charge::Number
    hist::Vector{Tuple{Number,Number}}
    annihilation_time::Union{Float64,Nothing}
    creation_time::Float64
    id_annihilator::Union{Int,Nothing}

    Defect(;id,charge,loc,t) = new(id,charge,[loc],nothing,t,nothing)
end
last_loc(d::Defect) = d.hist[end]
creation_loc(d::Defect) = d.hist[1]
push_position!(d::Defect,loc) = push!(d.hist,loc)

mutable struct DefectTracker
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsN::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Float64 # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(thetas,model,lattice) # constructor
        vortices,antivortices = spot_defects(thetas,model,lattice)
        defectsP = [Defect(id=i,charge=vortices[i][3],loc=vortices[i][1:2],t=model.t) for i in each(vortices)]
        defectsN = [Defect(id=i,charge=antivortices[i][3],loc=antivortices[i][1:2],t=model.t) for i in each(antivortices)]
        new(defectsP,defectsN,model.t)
    end
end
number_defectsP(dt::DefectTracker) = length(dt.defectsP)
number_defectsN(dt::DefectTracker) = length(dt.defectsN)
number_defects(dt::DefectTracker)  = length(dt.defectsN)  + length(dt.defectsN)
number_active_defectsP(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsP])
number_active_defectsN(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsN])
number_active_defects(dt::DefectTracker)  = number_active_defectsN(dt) + number_active_defectsP(dt)

function t_bounds(dft::DefectTracker)
    alldefects = vcat(dft.defectsP,dft.defectsN)
    tmin = minimum([d.creation_time for d in alldefects])
    tmax = maximum([d.annihilation_time for d in alldefects])
    return tmin, tmax
end

function ID_active_defects(dt::DefectTracker)
    activeP = Int[]
    for i in 1:number_defectsP(dt)
        if dt.defectsP[i].annihilation_time == nothing push!(activeP,i) end
    end
    activeN = Int[]
    for i in 1:number_defectsN(dt)
        if dt.defectsN[i].annihilation_time == nothing push!(activeN,i) end
    end
    return activeP,activeN
end

function add_defect!(dt::DefectTracker;charge,loc)
    if charge > 0 push!(dt.defectsP,Defect(id=1+number_defectsP(dt),charge=charge,loc=loc,t=dt.current_time))
    else          push!(dt.defectsN,Defect(id=1+number_defectsN(dt),charge=charge,loc=loc,t=dt.current_time))
    end
end

function pair_up_hungarian(dt::DefectTracker,new,old,lattice::AbstractLattice,charge::String)
    # charge can be either "+" or "-"
    distance_matrixx = distance_matrix(new,old,lattice) # m_new lignes, m_old colonnes
    proposal         = hungarian(distance_matrixx)[1] # length = length(new)
    assignment       = copy(proposal) # because it will be modified in the next for loop

    #= A few comments :
    1. The proposal is the match between indices of the vectors
    new,old while the assignment matches actual IDs of the DefectTracker.
    2. If cost_matrix is a NxM matrix (workers x jobs), the output of hungarian(cost_matrix)
    is a Nx1 vector containing the assignments of each of the N workers to the indices of the jobs (1 < indices < M).
    A "0" means the worker has not been assigned to any job, hence the assignment[i] ≠ 0 condition below.
    =#
    if charge == "+"
        for i in eachindex(assignment)
            for j in 1:number_defectsP(dt)
                if assignment[i] ≠ 0 && dt.defectsP[j].annihilation_time == nothing && last_loc(dt.defectsP[j]) == old[proposal[i]]
                    # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
                    # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
                    assignment[i] = j
                    break # breaks innerloop only
                end
            end
        end
    else # if charge  == "-"
        for i in eachindex(assignment)
            for j in 1:number_defectsN(dt)
                if assignment[i] ≠ 0 && dt.defectsN[j].annihilation_time == nothing && last_loc(dt.defectsN[j]) == old[proposal[i]]
                    assignment[i] = j
                    break
                end
            end
        end
    end
    return assignment
end

function find_closest_before_annihilation(dt,lattice,old_loc_defect)
    distance = Inf ; ID_antidefect = -1 # dummy
    for i in each(dt.defectsN)
        if dt.defectsN[i].annihilation_time == dt.current_time # it has just annihilated
            tmp = dist(lattice,old_loc_defect,last_loc(dt.defectsN[i]))
            if tmp < distance
                distance = tmp
                ID_antidefect = i
            end
        end
    end
    if ID_antidefect == -1
        println("Something weird is going on...")
        println([dt.defectsN[i].annihilation_time == dt.current_time for i in 1:number_defectsN(dt)])
    end
    return ID_antidefect,last_loc(dt.defectsN[ID_antidefect])
end

# estimation_location_annihilation() has been replaced by mean_2_positions()
function annihilate_defects(dt::DefectTracker,ids_annihilated_defects,L)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex,old_loc_antivortex = find_closest_before_annihilation(dt,lattice,old_loc_vortex)

        dt.defectsP[i].id_annihilator = ID_antivortex
        dt.defectsN[ID_antivortex].id_annihilator = i

        # dt.defectsP[i].annihilation_time = dt.current_time
        # dt.defectsN[ID_antivortex].annihilation_time = dt.current_time

        estimate = mean_2_positions(old_loc_vortex,old_loc_antivortex,L)
        push_position!(dt.defectsP[i],estimate)
        push_position!(dt.defectsN[ID_antivortex],estimate)
    end
    return dt
end


function update_and_track!(thetas::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice,dft::DefectTracker,tmax::Number,every::Number) where T<:AbstractFloat
    next_tracking_time = model.t
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ next_tracking_time || model.t ≥ tmax # second condition to end at the same time than the model
            next_tracking_time += every
            dft.current_time = model.t
            update_DefectTracker!(dft,thetas,model,lattice)
        end
    end
    return dft
end

function update_and_track_plot!(thetas::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice,dft::DefectTracker,tmax::Number,every::Number;defects=false) where T<:AbstractFloat
    next_tracking_time = model.t
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ next_tracking_time || model.t ≥ tmax # second condition to end at the same time than the model
            next_tracking_time += every
            dft.current_time = model.t
            update_DefectTracker!(dft,thetas,model,lattice)
            display(plot_thetas(thetas,model,lattice,defects=true))
        end
    end
    return dft
end

function update_DefectTracker!(dt::DefectTracker,thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::AbstractLattice)

    NB = lattice.periodic
    vortices_new,antivortices_new = spot_defects(thetas,model,lattice)

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    locP_old    = [last_loc(dt.defectsP[i]) for i in each(dt.defectsP)]
    locN_old    = [last_loc(dt.defectsN[i]) for i in each(dt.defectsN)]
    chargeP_old = [dt.defectsP[i].charge    for i in each(dt.defectsP)]
    chargeN_old = [dt.defectsN[i].charge    for i in each(dt.defectsN)]

    locP_new    = [vortices_new[i][1:2]     for i in each(vortices_new)]
    locN_new    = [antivortices_new[i][1:2] for i in each(antivortices_new)]
    chargeP_new = [vortices_new[i][3]       for i in each(vortices_new)]
    chargeN_new = [antivortices_new[i][3]   for i in each(antivortices_new)]

    Np_new,Np_old = length(locP_new),length(locP_old)
    Nn_new,Nn_old = length(locN_new),length(locN_old)
    N_old = Np_old + Nn_old
    N_new = Np_new + Nn_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nn_new == Nn_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        for i in 1:Np_new push_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end

    elseif Np_new == Np_old == 0 && Nn_new == Nn_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        for i in 1:Nn_new push_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new add_defect!(dt,charge=chargeP_new[i],loc=locP_new[i]) end
        for i in 1:Nn_new add_defect!(dt,charge=chargeN_new[i],loc=locN_new[i]) end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP,id_just_annihilated_defectM = ID_active_defects(dt) # seek for not yet annihilated defects
        # println("hello")
        # println(id_just_annihilated_defectP)
        # println(id_just_annihilated_defectM)
        for i in id_just_annihilated_defectP dt.defectsP[i].annihilation_time = dt.current_time end
        for i in id_just_annihilated_defectM dt.defectsN[i].annihilation_time = dt.current_time end
        # annihilation_time is already taken care of in the annihilate_defects function
        dt = annihilate_defects(dt::DefectTracker,id_just_annihilated_defectP,L)

    elseif Np_new > 0 && Np_old > 0 && Nn_old > 0 && Nn_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices) push_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsP(dt)
            if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices,i)
            end
        end
        ID_annihilated_antivortices = ID_active_defects(dt)[2] # in this special case, there is no antivortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end
        dt = annihilate_defects(dt,ID_annihilated_vortices,L)

    elseif Nn_new > 0 && Nn_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices) push_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:number_defectsN(dt)
            if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_antivortices,i)
            end
        end
        ID_annihilated_vortices = ID_active_defects(dt)[1] # in this special case, there is no vortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = dt.current_time end
        for i in ID_annihilated_antivortices dt.defectsN[i].annihilation_time = dt.current_time end

        dt = annihilate_defects(dt,ID_annihilated_vortices,L)
    else # end of special cases

    # GENERAL TREATMENT
        assignment_vortices     = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new push_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
            for i in 1:Nn_new push_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

        # CASE 2 : creation !
    elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in each(loc_created_vortex) add_defect!(dt,charge=chargeP_new[j],loc=loc_created_vortex[j][1:2]) end

            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in each(loc_created_antivortex) add_defect!(dt,charge=chargeN_new[j],loc=loc_created_antivortex[j][1:2]) end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    push_position!(dt.defectsP[assignment_vortices[i]],locP_new[i])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    push_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i])
                end
            end

        # CASE 3 : annihilation !
    elseif N_new < N_old
             # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
             for i in eachindex(assignment_vortices)     push_position!(dt.defectsP[assignment_vortices[i]],locP_new[i]) end
             for i in eachindex(assignment_antivortices) push_position!(dt.defectsN[assignment_antivortices[i]],locN_new[i]) end

            # Identify annihilated defects
            ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
            for i in 1:number_defectsP(dt)
                if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    dt.defectsP[i].annihilation_time = dt.current_time # from now on, defects that have just annihilated have annihilation_time == t
                    push!(ID_annihilated_vortices,i)
                end
            end
            for i in 1:number_defectsN(dt)
                if i ∉ assignment_antivortices && dt.defectsN[i].annihilation_time == nothing
                    dt.defectsN[i].annihilation_time = dt.current_time
                    push!(ID_annihilated_antivortices,i)
                end
            end
            # @assert length(ID_annihilated_vortices) == length(ID_annihilated_antivortices)
            # println(length(assignment_vortices)," ",length(assignment_antivortices))
            dt = annihilate_defects(dt,ID_annihilated_vortices,L)
        end # end of general treatment
    end # end of special cases & general treatment
    return dt
end

function MSD(dft::DefectTracker,model::AbstractModel,lattice::AbstractLattice)
    nP = number_defectsP(dft)
    nN = number_defectsN(dft)
    tmin,tmax = t_bounds(dft) # (tmin,tmax) = timestamps of (first defect creation , last defect annihilation)

    hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1

    # Compute the SD
    SD_P = NaN*zeros(nP,Int(tmax))
    SD_N = NaN*zeros(nN,Int(tmax))
    for n in 1:nP
        defect = dft.defectsP[n]
        index_creation = round(Int,defect.creation_time/dummy_dt)
        index_annihilation = round(Int,defect.annihilation_time/dummy_dt)
        SD_P[n,index_creation:index_annihilation] = square_displacement(defect,lattice)
    end
    for n in 1:nN
        defect = dft.defectsN[n]
        index_creation = round(Int,defect.creation_time/dummy_dt)
        index_annihilation = round(Int,defect.annihilation_time/dummy_dt)
        SD_N[n,index_creation:index_annihilation] = square_displacement(defect,lattice)
    end

    # Now average to compute the MSD
    MSD_P = nanmean(SD_P,1)[1,:]
    MSD_N = nanmean(SD_N,1)[1,:]
    MSD_N = nanmean(hcat(SD_M,SD_N),1)[1,:]

    return MSD, MSD_P, MSD_N
end

function square_displacement(d::Defect,lattice::AbstractLattice)
    loc_t0 = creation_loc(d)
    return [dist(lattice,pos,loc_t0) for loc in d.hist] .^ 2
end

function interdefect_distance(dft,defect1,defect2,lattice)
    # TODO take care of case with creation and/or annihilation time different.
    # So far, this care is left to the user...
    @assert defect1.creation_time == defect2.creation_time
    @assert defect1.annihilation_time == defect2.annihilation_time

    R = [dist(lattice,defect1.hist[t],defect2.hist[t]) for t in each(defect1.hist)]
    return R
end
## Small helpful methods for scripts
function number_defects(model::AbstractModel,lattice::AbstractLattice)
    a,b = spot_defects(thetas,T,BC)
    return length(a) + length(b)
end

function theta_mid(x::T,y::T,symm)::T  where T<:AbstractFloat
    arcl = arclength(x,y,symm)
    return x + arcl/2 , y - arcl/2
end
