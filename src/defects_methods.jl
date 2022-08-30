include("../src/lattices.jl");
include("../src/models.jl");
include("../src/core_methods.jl");

function arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat
    #= This function returns the signed arclength on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
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

function get_vorticity(thetasmodpi::Matrix{T},model::AbstractModel{T},lattice::AbstractLattice,i::Int,j::Int)::T where T<:AbstractFloat
    #= Note : thetasmodpi = mod.(thetas,symm)
        By convention, the neighbours have the following ordering
           3   2
        4    x   1
           5   6
        =#
    symm = sym(model)
    angles_corners = get_neighbours(thetasmodpi,model,lattice,i,j,is_in_bulk(i,j,lattice.L))
    perimeter_covered = 0.0
    for i in 1:length(angles_corners)-1
        perimeter_covered += arclength(angles_corners[i],angles_corners[i+1],symm)
    end
    perimeter_covered += arclength(angles_corners[end],angles_corners[1],symm)
    charge = round(perimeter_covered/2π,digits=1) # here 2π no matter the symmetry of the model. if nematic, will return half-integer defects
    return charge
end

function spot_defects(thetas,model::AbstractModel,lattice::AbstractLattice)
    L = model.L
    if lattice.periodic range_bc = 1:L else range_bc = 2:L-1 end
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

    seuil = 4

    elements_to_keep_minus = trues(length(list_vortices_minus))
    for i in 2:length(list_vortices_minus)
        for j in 1:i-1
            if dist(lattice,list_vortices_minus[i],list_vortices_minus[j]) ≤ seuil
                elements_to_keep_minus[i] = false
                # break # will only break the inner for loop and continue the outer one
            end
        end
    end
    elements_to_keep_plus  = trues(length(list_vortices_plus))
    for i in 2:length(list_vortices_plus)
        for j in 1:i-1
            if dist(lattice,list_vortices_plus[i],list_vortices_plus[j]) ≤ seuil
                elements_to_keep_plus[i] = false
                # break # will only break the inner for loop and continue the outer one
            end
        end
    end
    # return list_vortices_plus,list_vortices_minus
    return list_vortices_plus[elements_to_keep_plus],list_vortices_minus[elements_to_keep_minus]
end

#is it even useful now ?
# function spot_single_default_global(thetas::Matrix{T})::Tuple{Tuple{Int16,Int16},T} where T<:AbstractFloat
#     L = size(thetas)[1]
#     list_defaults = Tuple{Int16,Int16,T}[]
#     for i in 2:L-1
#         for j in 2:L-1
#             q = get_vorticity(thetas,i,j,L)
#             if abs(q) > 0.1
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
    charge::Float32
    hist::Vector{Tuple{Float32,Float32}}
    annihilation_time::Union{Nothing,Int}
    creation_time::Int
    id_annihilator::Union{Nothing,Int}

    function Defect(;id,charge,loc,t) # constructor
        new(id,charge,[loc],nothing,t,nothing)
    end
end
    function last_loc(d::Defect)     return d.hist[end] end
    function creation_loc(d::Defect) return d.hist[1] end

mutable struct DefectTracker
    Np::Int # length = number of (+) defects
    Nm::Int # length = number of (-) defects
    Npactive::Int
    Nmactive::Int
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsM::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Int # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(;thetas,T,BC,t) # constructor
        vortices,antivortices = spot_defects(thetas,T,BC)
        Np = length(vortices) ; Nm = length(antivortices)
        new(Np,Nm,Np,Nm,[Defect(id=i,charge=+1/2,loc=vortices[i],t=t) for i in 1:Np],[Defect(id=i,charge=-1/2,loc=antivortices[i],t=t) for i in 1:Nm],t)
    end
end

function total_tracked(dt::DefectTracker) return dt.Np + dt.Nm end

function SD(dt::DefectTracker,L)
        lengthsP  = [length(dt.defectsP[n].hist) for n in 1:dt.Np]
        lengthsM  = [length(dt.defectsM[n].hist) for n in 1:dt.Nm]
        maxlength = maximum(vcat(lengthsP,lengthsM))
        SDp = NaN*zeros(dt.Np,maxlength) # square displacement for positive defects
        SDm = NaN*zeros(dt.Nm,maxlength) # square displacement for negative defects

        for n in 1:dt.Np
            SDp[n,1:lengthsP[n]] = ([dist(dt.defectsP[n].hist[i] , creation_loc(dt.defectsP[n]),L,squared=true) for i in 1:lengthsP[n]])
            # SDp[n,1:lengthsP[n]] = reverse([dist(dt.defectsP[n].hist[i] , last_loc(dt.defectsP[n]),L,squared=true) for i in 1:lengthsP[n]])
        end
        for n in 1:dt.Nm
            SDm[n,1:lengthsM[n]] = ([dist(dt.defectsM[n].hist[i] , creation_loc(dt.defectsM[n]),L,squared=true) for i in 1:lengthsM[n]])
            # SDm[n,1:lengthsM[n]] = reverse([dist(dt.defectsM[n].hist[i] , last_loc(dt.defectsM[n]),L,squared=true) for i in 1:lengthsM[n]])
        end


        MSDp = nanmean(SDp,1)[1,:]
        MSDm = nanmean(SDm,1)[1,:]
        MSD  = (MSDm .+ MSDp)/2

        return SDp,SDm,MSDp,MSDm,MSD
    end

function add_defect!(dt::DefectTracker;charge,loc,t)
    if charge > 0
        dt.Np += 1
        push!(dt.defectsP,Defect(id=dt.Np,charge=charge,loc=loc,t=t))
    else
        dt.Nm += 1
        push!(dt.defectsM,Defect(id=dt.Nm,charge=charge,loc=loc,t=t))
    end
end

# TODO le code pourrait etre bien plus simple, avec un simple count + condition
function defects_active(dt::DefectTracker)
    P = [] ; M = []
    for n in 1:dt.Np
        if dt.defectsP[n].annihilation_time == nothing push!(P,n) end
    end
    for n in 1:dt.Nm
        if dt.defectsM[n].annihilation_time == nothing push!(M,n) end
    end
    return P,M,length(P)+length(M)
end

function pair_up_hungarian(dt::DefectTracker,new,old,L,charge)
    distance_matrixx    = distance_matrix(new,old,L) # m_new lignes, m_old colonnes
    proposal            = hungarian(distance_matrixx)[1] # length = length(new)
    assignment          = copy(proposal) # because it will be modified in the next for loop

    #= A few comments :
    1. The proposal is the match between indices of the vectors
    new,old while the assignment matches actual IDs of the DefectTracker.
    2. If cost_matrix is a NxM matrix (workers x jobs), the output of hungarian(cost_matrix)
    is a Nx1 vector containing the assignments of each of the N workers to the indices of the jobs (1 < indices < M).
    A "0" means the worker has not been assigned to any job, hence the assignment[i] ≠ 0 condition below.
    =#
    if charge > 0
        for i in eachindex(assignment)
            for j in 1:dt.Np
                if assignment[i] ≠ 0 && dt.defectsP[j].annihilation_time == nothing && last_loc(dt.defectsP[j]) == old[proposal[i]]
                    # the first condition is a protection against the creation case, where the Hungarian algo matches the newly created vortex to 0
                    # the second condition ensures that one only considers currently living vortices and not vortices now annihilated
                    assignment[i] = j
                    break # breaks innerloop only
                end
            end
        end
    else # if charge < 0
        for i in eachindex(assignment)
            for j in 1:dt.Nm
                if assignment[i] ≠ 0 && dt.defectsM[j].annihilation_time == nothing && last_loc(dt.defectsM[j]) == old[proposal[i]]
                    assignment[i] = j
                    break
                end
            end
        end
    end
    return assignment
end

function find_closest_before_annihilation(dt,old_loc_defect,L)
    distance = Inf ; ID_antidefect = -1 # dummy
    for i in each(dt.defectsM)
        if dt.defectsM[i].annihilation_time == dt.current_time # it has just annihilated
            tmp = dist(old_loc_defect,last_loc(dt.defectsM[i]),L)
            if tmp < distance
                distance = tmp
                ID_antidefect = i
            end
        end
    end
    if ID_antidefect == -1
        println([dt.defectsM[i].annihilation_time == dt.current_time for i in 1:dt.Nm])
    end
    return ID_antidefect,last_loc(dt.defectsM[ID_antidefect])
end

function estimation_location_annihilation((a,b),(x,y),L)
    dx = (x - a) #; dx = min(dx,L-dx)
    if abs(L-dx) < abs(dx) dx = -(L-dx) end
    # if L-dx < dx dx = -(L-dx) end
    dy = (y - b) #; dy = min(dy,L-dy)
    if abs(L-dy) < abs(dy) dy = -(L-dy) end
    # if L-dy < dy dy = -(L-dy) end
    estimate_loc_collision = mod1.((a,b) .+ 0.5.*(dx,dy),L)
    # estimate_loc_collision = mod1.(0.5.*(a,b) .+ 0.5.*(x,y),L)
    return estimate_loc_collision
end
# estimation_location_annihilation((50,50),(60,60),L) == (55,55)
# estimation_location_annihilation((10,10),(90,90),L) == (100,100)
# estimation_location_annihilation((49,66),(51,61),L) == (50.0, 63.5)

function annihilate_defects(dt::DefectTracker,ids_annihilated_defects,L)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex,old_loc_antivortex = find_closest_before_annihilation(dt,old_loc_vortex,L)

        dt.defectsP[i].id_annihilator = ID_antivortex
        dt.defectsM[ID_antivortex].id_annihilator = i

        # dt.defectsP[i].annihilation_time = dt.current_time
        # dt.defectsM[ID_antivortex].annihilation_time = dt.current_time

        estimate = estimation_location_annihilation(old_loc_vortex,old_loc_antivortex,L)
        push!(dt.defectsP[i].hist,estimate)
        push!(dt.defectsM[ID_antivortex].hist,estimate)
    end
    return dt
end

function update_DefectTracker(dt::DefectTracker,thetas::Matrix{Float32},BC,
    vortices_new::Vector{Tuple{Int,Int}},antivortices_new::Vector{Tuple{Int,Int}},
    vortices_old::Vector{Tuple{Int,Int}},antivortices_old::Vector{Tuple{Int,Int}},t)

    dt.current_time = t

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    Np_new,Np_old = length(vortices_new),length(vortices_old)
    Nm_new,Nm_old = length(antivortices_new),length(antivortices_old)
    N_old = Np_old + Nm_old ; N_new = Np_new + Nm_new

    dt.Npactive = Np_new
    dt.Nmactive = Nm_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nm_new == Nm_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        assignment_vortices     = pair_up_hungarian(dt,vortices_new,vortices_old,L,+1/2)
        for i in 1:Np_new push!(dt.defectsP[assignment_vortices[i]].hist,vortices_new[i]) end

    elseif Np_new == Np_old == 0 && Nm_new == Nm_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dt,antivortices_new,antivortices_old,L,-1/2)
        for i in 1:Nm_new push!(dt.defectsM[assignment_antivortices[i]].hist,antivortices_new[i]) end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new add_defect(dt,charge=+1/2,loc=vortices_new[i],t=t) end
        for i in 1:Nm_new add_defect(dt,charge=-1/2,loc=antivortices_new[i],t=t) end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP,id_just_annihilated_defectM,~ = defects_alive(dt) # seek for not yet annihilated defects
        for i in id_just_annihilated_defectP dt.defectsP[i].annihilation_time = t end
        for i in id_just_annihilated_defectM dt.defectsM[i].annihilation_time = t end
        # annihilation_time is already taken care of in the annihilate_defects function
        dt = annihilate_defects(dt::DefectTracker,id_just_annihilated_defectP,L)

    elseif Np_new > 0 && Np_old > 0 && Nm_old > 0 && Nm_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dt,vortices_new,vortices_old,L,+1/2)
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices) push!(dt.defectsP[assignment_vortices[i]].hist,vortices_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:dt.Np
            if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_vortices,i)
            end
        end
        ID_annihilated_antivortices = defects_alive(dt)[2] # in this special case, there is no antivortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = t end
        for i in ID_annihilated_antivortices dt.defectsM[i].annihilation_time = t end
        dt = annihilate_defects(dt,ID_annihilated_vortices,L)

    elseif Nm_new > 0 && Nm_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,antivortices_new,antivortices_old,L,-1/2)
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices) push!(dt.defectsM[assignment_antivortices[i]].hist,antivortices_new[i]) end
        # Identify annihilated defects
        ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
        for i in 1:dt.Nm
            if i ∉ assignment_antivortices && dt.defectsM[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                push!(ID_annihilated_antivortices,i)
            end
        end
        ID_annihilated_vortices = defects_alive(dt)[1] # in this special case, there is no vortices left in the "new" timestep

        for i in ID_annihilated_vortices     dt.defectsP[i].annihilation_time = t end
        for i in ID_annihilated_antivortices dt.defectsM[i].annihilation_time = t end

        dt = annihilate_defects(dt,ID_annihilated_vortices,L)
    else # end of special cases
    # GENERAL TREATMENT
        assignment_vortices     = pair_up_hungarian(dt,vortices_new,vortices_old,L,+1/2)
        assignment_antivortices = pair_up_hungarian(dt,antivortices_new,antivortices_old,L,-1/2)

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new push!(dt.defectsP[assignment_vortices[i]].hist,vortices_new[i]) end
            for i in 1:Nm_new push!(dt.defectsM[assignment_antivortices[i]].hist,antivortices_new[i]) end

        # CASE 2 : creation !
    elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in each(loc_created_vortex) add_defect(dt,charge=+1/2,loc=loc_created_vortex[j],t=t) end

            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in each(loc_created_antivortex) add_defect(dt,charge=-1/2,loc=loc_created_antivortex[j],t=t) end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    push!(dt.defectsP[assignment_vortices[i]].hist,vortices_new[i])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    push!(dt.defectsM[assignment_antivortices[i]].hist,antivortices_new[i])
                end
            end

        # CASE 3 : annihilation !
    elseif N_new < N_old
             # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
             for i in eachindex(assignment_vortices)     push!(dt.defectsP[assignment_vortices[i]].hist,vortices_new[i]) end
             for i in eachindex(assignment_antivortices) push!(dt.defectsM[assignment_antivortices[i]].hist,antivortices_new[i]) end

            # Identify annihilated defects
            ID_annihilated_vortices = [] ; ID_annihilated_antivortices = []
            for i in 1:dt.Np
                if i ∉ assignment_vortices && dt.defectsP[i].annihilation_time == nothing # a vortex has just annihilated if its ID is not in the assignment list AND if its annihilation time is still "nothing"
                    dt.defectsP[i].annihilation_time = t # from now on, defects that have just annihilated have annihilation_time == t
                    push!(ID_annihilated_vortices,i)
                end
            end
            for i in 1:dt.Nm
                if i ∉ assignment_antivortices && dt.defectsM[i].annihilation_time == nothing
                    dt.defectsM[i].annihilation_time = t
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

## Small helpful methods for scripts

function number_defects(model::AbstractModel,lattice::AbstractLattice)
    a,b = spot_defects(thetas,T,BC)
    return length(a) + length(b)
end

function theta_mid(x::T,y::T,symm)::T  where T<:AbstractFloat
    arcl = arclength(x,y,symm)
    return x + arcl/2 , y - arcl/2
end
