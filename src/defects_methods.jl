include("lattices.jl");
include("models.jl");
include("core_methods.jl");

using Hungarian
import Random.shuffle

function arclength(theta1::T,theta2::T,symm)::T where T<:AbstractFloat
    #= This function returns the signed arclength (in radians) on the unit trigonometric circle .
    Clockwise        >> sign -
    Counterclockwise >> sign +
    WARNING : Note that the inputs thetas need to lie within [0,π] or [0,2π], depending on the symmetry "symm" of the model =#
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

function get_vorticity(thetasmod::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice,i::Int,j::Int)::T where T<:AbstractFloat
    #= Note : thetasmod = mod.(thetas,symm*pi)
        By convention, for a SquareLattice, the neighbours have the following ordering
           2                            3
        3  x  1 when plotted and     4  x  2  in the matrix form
           4                            1
       By convention, for a TriangularLattice, the neighbours have the following ordering
          3  2                         4  3
        4  x  1 when plotted and     5  x  2  in the matrix form
          5  6                         6  1

        Indeed, recall that the visualisation procedure operates a 90° counterclockwise rotation.
        We thus anticipate the rotation already in the matrix form, so that the first neighbour
        is on the right, following mathematical use (associated to angle θ = 0).
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

number_defects(thetas,model,lattice) = sum(length,spot_defects(thetas,model,lattice))
function spot_defects(thetas::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice;find_type=false) where T<:AbstractFloat
    L = lattice.L
    vortices_plus  = Tuple{Int,Int,T,T}[]
    vortices_minus = Tuple{Int,Int,T,T}[]

    if     model.symmetry == "nematic" modd = T(pi)
    elseif model.symmetry == "polar"   modd = T(2pi)
    end
    thetasmod = mod.(copy(thetas),modd)
    if model.rho < 1 preconditionning!(thetasmod,model,lattice) end
    # trelax = 0.3 ; relax!(thetasmod,model,trelax)
    if lattice.periodic range_bc = 1:L else range_bc = 2:L-1 end
    temporary_type = NaN
    for i in range_bc
        for j in range_bc
            q = get_vorticity(thetasmod,model,lattice,i,j)
            if     q > + 0.1 push!(vortices_plus,(i,j,q,temporary_type)) # we want to keep ±½ defects, and not rounding errors
            elseif q < - 0.1 push!(vortices_minus,(i,j,q,temporary_type))
            end
        end
    end
    vortices_plus_no_duplicates  = merge_duplicates(vortices_plus,lattice)
    vortices_minus_no_duplicates = merge_duplicates(vortices_minus,lattice)

    if find_type return find_types(vortices_plus_no_duplicates,vortices_minus_no_duplicates,mod.(thetas,2pi),lattice)
    else return vortices_plus_no_duplicates,vortices_minus_no_duplicates
    end
end

function find_types(list_p,list_n,thetas,lattice)
    # Positive defects
    pos_p    = [list_p[i][1:2] for i in each(list_p)]
    charge_p = [list_p[i][3]   for i in each(list_p)]
    type_p   = [list_p[i][4]   for i in each(list_p)]
    # Negative defects
    pos_n    = [list_n[i][1:2] for i in each(list_n)]
    charge_n = [list_n[i][3]   for i in each(list_n)]
    type_n   = [list_n[i][4]   for i in each(list_n)]
    # All defects together
    pos_all = vcat(pos_n,pos_p)
    type_all = vcat(type_n,type_p)

    # define Denoising AutoEncoder DAE
    if length(charge_p) > 0
        if charge_p[1] == 1
            DAE = DAE_positive1
        elseif charge_p[1] == 0.5
            DAE = DAE_positive12
        end
    end

    total_number_defects = length(pos_n) + length(pos_p)
    density_defects = total_number_defects / lattice.L^2

    if true #density_defects < 1/(2WINDOW+1)^2
        #= WINDOW x WINDOW portion around the defect and we want this square
        free of any other defect to proceed. If too crowded, don't
        even enter the computationally expensive operations hereafter.
        =#
        for n in each(pos_p)
            if alone_in_window(pos_p[n],pos_all,lattice,WINDOW) # heavy, computes distance
                i,j = pos_p[n]
                no_problem_go_ahead,thetas_zoom = zoom(thetas,lattice,i,j,WINDOW)
                if no_problem_go_ahead
                    #= A problem could occur if defect close to boundary
                    and lattice not periodic. If so, leave the type value
                    unchanged, i.e =NaN =#
                    if charge_p[n] == +1
                        denoised_theta_zoom = DAE(reshape(provide_div_rot(thetas_zoom),(W21,W21,3,1)))
                        denoised_theta_zoom_reshaped = reshape(denoised_theta_zoom,(W21,W21))
                        type_p[n] = infer_mu(denoised_theta_zoom_reshaped,q=charge_p[n])
                    end
                end
            end
        end
        for n in each(pos_n)
            if alone_in_window(pos_n[n],pos_all,lattice,WINDOW) # heavy, computes distance
                i,j = pos_n[n]
                no_problem_go_ahead,thetas_zoom = zoom(thetas,lattice,i,j,WINDOW)
                if no_problem_go_ahead
                    #= A problem could occur if defect close to boundary
                    and lattice not periodic. If so, leave the type value
                    unchanged, i.e "unknown" =#
                    if charge_n[n] == -1
                        type_n[n] = infer_mu(thetas_zoom,q=charge_n[n])
                    end
                end
            end
        end
    end

    list_p_updated = similar(list_p)
    for i in each(list_p) list_p_updated[i] = (list_p[i][1],list_p[i][2],list_p[i][3],type_p[i]) end
    list_n_updated = similar(list_n)
    for i in each(list_n) list_n_updated[i] = (list_n[i][1],list_n[i][2],list_n[i][3],type_n[i]) end

    return list_p_updated,list_n_updated
end

function alone_in_window(pos,pos_all,lattice,window)::Bool
    alone_in_window = true
    for loc in pos_all
        distance = dist(lattice,pos,loc)
        if 0 < distance ≤ window+1 # take margin because non integer positions might get rounded
            alone_in_window = false
            break
        end
    end
    return alone_in_window
end

function merge_duplicates(list,lattice;radius=4)
    #= In this list, there might be doubles/triples (2/3 locations for the
    same physical vortex). We thus seek for numerically identified vortices
    which are neighbours and with the same charge to delete them. =#
    pos    = [list[i][1:2] for i in each(list)]
    charge = [list[i][3]   for i in each(list)]
    type   = [list[i][4]   for i in each(list)]
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
            mean_loc_defect = mean_N_positions(tmp,lattice.L,true)
            push!(merged_duplicates,(mean_loc_defect[1],mean_loc_defect[2],charge[i],type[i]))
        end
    end
    return merged_duplicates
end

function preconditionning!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Abstract2DLattice)
    remove_isolates!(thetas,model,lattice)
    fill_holes!(thetas,model,lattice)
    trelax = 0.5
    relax!(thetas,model,trelax)
end

function remove_isolates!(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Abstract2DLattice)
    L = size(thetas,1)
    for j in 1:L
        for i in 1:L
            if isempty(get_neighbours(thetas,model,lattice,i,j,is_in_bulk(i,j,L)))
                thetas[i,j] = NaN
            end
        end
    end
    return thetas
end

function fill_holes!(thetas::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice) where T<:AbstractFloat
    L = size(thetas,1)
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

# function vortices_to_keep(vortices,lattice,seuil)
#     elements_to_keep = trues(length(vortices))
#     for i in 2:length(vortices)
#         for j in 1:i-1
#             if dist(lattice,vortices[i][1:2],vortices[j][1:2]) ≤ seuil
#                 elements_to_keep[i] = false
#                 # break # will only break the inner for loop and continue the outer one
#             end
#         end
#     end
#     return elements_to_keep
# end

## Defect Tracking (with a lot of defects, creation and annihilation)
mutable struct Defect
    id::Int
    charge::Number
    thetas_zoom::Vector{Matrix{Float32}} # types might change over the simulation
    type::Vector{Float64} # types might change over the simulation
    pos::Vector{Tuple{Number,Number}}
    annihilation_time::Union{Float64,Nothing}
    creation_time::Float64
    id_annihilator::Union{Int,Nothing}
end
Defect(;id,charge,loc,thetas_zoom,t,type=NaN) = Defect(id,charge,[thetas_zoom],[type],[loc],nothing,t,nothing)
# Defect(;id,charge,loc,t,type=NaN) = Defect(id,charge,[type],[loc],nothing,t,nothing)

first_loc(d::Defect) = d.pos[1]
last_loc(d::Defect)  = d.pos[end]

first_type(d::Defect) = d.type[1]
last_type(d::Defect)  = d.type[end]

# function update_position_and_type!(d::Defect,new_loc,new_type)
#     push!(d.pos,new_loc)
#     if new_type == "unknown"  push!(d.type,last_type(d)) # by default, if unknown, push the last known type
#     else push!(d.type,new_type)
#     end)
# end
function update_position_and_type!(d::Defect,new_loc,new_type,thetas_zoomed)
    push!(d.pos,new_loc)
    push!(d.type,new_type)
    push!(d.thetas_zoom,thetas_zoomed)
end

# n'a plus lieu d'être
# function number_type_changes(d::Defect)
#     number_changes = 0
#     for i in 2:length(d.type) if d.type[i] ≠ d.type[i-1] number_changes += 1 end end
#     return number_changes
# end

mutable struct DefectTracker
    defectsP::Vector{Defect} # the id of a defect is its index in this vector
    defectsN::Vector{Defect} # so there is a (+)defect with id=1 AND and a (-)defect with id=1
    current_time::Float64 # latest update time (by convention, the creation time of the whole data structure = 0)

    function DefectTracker(thetas,model,lattice;find_type=false) # constructor
        vortices,antivortices = spot_defects(thetas,model,lattice;find_type=find_type)
        if find_type
            defectsP = [Defect(id=i,charge=vortices[i][3],type=vortices[i][4],thetas_zoom=zoom(thetas,lattice,vortices[i][1:2]...,WINDOW)[2],loc=vortices[i][1:2],t=model.t) for i in each(vortices)]
            defectsN = [Defect(id=i,charge=antivortices[i][3],type=antivortices[i][4],thetas_zoom=zoom(thetas,lattice,antivortices[i][1:2]...,WINDOW)[2],loc=antivortices[i][1:2],t=model.t) for i in each(antivortices)]
        else
            defectsP = [Defect(id=i,charge=vortices[i][3],thetas_zoom=zoom(thetas,lattice,vortices[i][1:2]...,WINDOW)[2],type=NaN,loc=vortices[i][1:2],t=model.t) for i in each(vortices)]
            defectsN = [Defect(id=i,charge=antivortices[i][3],thetas_zoom=zoom(thetas,lattice,antivortices[i][1:2]...,WINDOW)[2],type=NaN,loc=antivortices[i][1:2],t=model.t) for i in each(antivortices)]
        end
        new(defectsP,defectsN,model.t)
    end
end
number_defectsP(dt::DefectTracker) = length(dt.defectsP)
number_defectsN(dt::DefectTracker) = length(dt.defectsN)
number_defects(dt::DefectTracker)  = length(dt.defectsN)  + length(dt.defectsN)
# function number_defects_type(dt::DefectTracker,type)
#     if     type in ["join","split","threefold1","threefold2"]
#         return count(isequal(type),[last_type(dt.defectsN[i]) for i in each(dt.defectsN)])
#     elseif type in ["source","sink","clockwise","counterclockwise"]
#         return count(isequal(type),[last_type(dt.defectsP[i]) for i in each(dt.defectsP)])
#     end
# end
# number_defects_types(dt::DefectTracker) = [number_defects_type(dt,type) for type in ["source","sink","clockwise","counterclockwise","join","split","threefold1","threefold2"]]
number_active_defectsP(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsP])
number_active_defectsN(dt::DefectTracker) = count(isnothing,[d.annihilation_time for d in dt.defectsN])
number_active_defects(dt::DefectTracker)  = number_active_defectsN(dt) + number_active_defectsP(dt)

first_types(dt::DefectTracker) = first_type.(vcat(dt.defectsP,dt.defectsN))
last_types(dt::DefectTracker) = last_type.(vcat(dt.defectsP,dt.defectsN))
types(dt::DefectTracker) = last_types(dt)

function t_bounds(dft::DefectTracker)
    alldefects = vcat(dft.defectsP,dft.defectsN)
    creation_times     = filter(!isnothing,[d.creation_time for d in alldefects])
    annihilation_times = filter(!isnothing,[d.annihilation_time for d in alldefects])
    if isempty(creation_times) tmin = 0.0
    else tmin = minimum(creation_times)
    end
    if isempty(annihilation_times) tmax = dft.current_time
    else tmax = maximum(creation_times)
    end
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

function add_defect!(dt::DefectTracker;charge,loc,type=NaN,thetas_zoom)
    if charge > 0 push!(dt.defectsP,Defect(id=1+number_defectsP(dt),charge=charge,thetas_zoom=thetas_zoom,type=type,loc=loc,t=dt.current_time))
    else          push!(dt.defectsN,Defect(id=1+number_defectsN(dt),charge=charge,thetas_zoom=thetas_zoom,type=type,loc=loc,t=dt.current_time))
    end
end

function pair_up_hungarian(dt::DefectTracker,new,old,lattice::Abstract2DLattice,charge::String)
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
        @warn "Annihilating defect not found: annihilation location will not be computed."
        # println("Something weird is going on...")
        println([dt.defectsN[i].annihilation_time == dt.current_time for i in 1:number_defectsN(dt)])
        return -1,(-1,-1) # dummy
    else
        return ID_antidefect,last_loc(dt.defectsN[ID_antidefect])
    end
end

function annihilate_defects(dt::DefectTracker,ids_annihilated_defects,lattice)
    for i in ids_annihilated_defects
        old_loc_vortex = last_loc(dt.defectsP[i])
        ID_antivortex,old_loc_antivortex = find_closest_before_annihilation(dt,lattice,old_loc_vortex)

        if ID_antivortex ≠ -1
            dt.defectsP[i].id_annihilator = ID_antivortex
            dt.defectsN[ID_antivortex].id_annihilator = i

            # dt.defectsP[i].annihilation_time = dt.current_time
            # dt.defectsN[ID_antivortex].annihilation_time = dt.current_time

            estimate = mean_2_positions(old_loc_vortex,old_loc_antivortex,lattice.L)
            update_position_and_type!(dt.defectsP[i],estimate,last_type(dt.defectsP[i]),NaN*zeros(WINDOW,WINDOW)) # dummy thetas_zoom full of NaN
            update_position_and_type!(dt.defectsN[ID_antivortex],estimate,last_type(dt.defectsN[i]),NaN*zeros(WINDOW,WINDOW)) # dummy thetas_zoom full of NaN
        end
    end
    return dt
end


function update_and_track!(thetas::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice,dft::DefectTracker,tmax::Number,every::Number;find_type=false,verbose=false) where T<:AbstractFloat
    next_tracking_time = model.t
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ next_tracking_time || model.t ≥ tmax # second condition to end at the same time than the model
            if verbose println("t = ",round(next_tracking_time,digits=1)) end
            next_tracking_time += every
            update_DefectTracker!(dft,thetas,model,lattice,find_type=find_type)
        end
    end
    return dft
end

function update_and_track_plot!(thetas::Matrix{T},model::AbstractModel{T},lattice::Abstract2DLattice,dft::DefectTracker,tmax::Number,every::Number;defects=false,find_type=false) where T<:AbstractFloat
    next_tracking_time = model.t
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ next_tracking_time || model.t ≥ tmax # second condition to end at the same time than the model
            println("t = ",model.t)
            next_tracking_time += every
            update_DefectTracker!(dft,thetas,model,lattice,find_type=find_type)
            display(plot_thetas(thetas,model,lattice,defects=defects))
        end
    end
    return dft
end

function update_DefectTracker!(dt::DefectTracker,thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::Abstract2DLattice;find_type=false)
    dt.current_time = model.t
    # NB = lattice.periodic
    vortices_new,antivortices_new = spot_defects(thetas,model,lattice,find_type=find_type)

    # if BC == "periodic" @assert length(vortices_new) == length(antivortices_new) && length(vortices_old) == length(antivortices_old) end
    locP_old    = [last_loc(dt.defectsP[i]) for i in each(dt.defectsP)]
    locN_old    = [last_loc(dt.defectsN[i]) for i in each(dt.defectsN)]
    chargeP_old = [dt.defectsP[i].charge    for i in each(dt.defectsP)]
    chargeN_old = [dt.defectsN[i].charge    for i in each(dt.defectsN)]

    locP_new    = [vortices_new[i][1:2]     for i in each(vortices_new)]
    locN_new    = [antivortices_new[i][1:2] for i in each(antivortices_new)]
    chargeP_new = [vortices_new[i][3]       for i in each(vortices_new)]
    chargeN_new = [antivortices_new[i][3]   for i in each(antivortices_new)]
    typeP_new   = [vortices_new[i][4]       for i in each(vortices_new)]
    typeN_new   = [antivortices_new[i][4]   for i in each(antivortices_new)]

    Np_new,Np_old = length(locP_new),length(locP_old)
    Nn_new,Nn_old = length(locN_new),length(locN_old)
    N_old = Np_old + Nn_old
    N_new = Np_new + Nn_new

    # Special simple cases to deal with upstream
    if N_new == N_old == 0 # do nothing, otherwise, "reducing over empty collection blablabla"

    elseif Nn_new == Nn_old == 0 && Np_new == Np_old > 0 # there are only (+) defects and no creation/annihilation
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        for i in 1:Np_new update_position_and_type!(dt.defectsP[assignment_vortices[i]],locP_new[i],typeP_new[i],zoom(thetas,lattice,locP_new[i]...,WINDOW)[2]) end

    elseif Np_new == Np_old == 0 && Nn_new == Nn_old > 0 # there are only (-) defects and no creation/annihilation
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        for i in 1:Nn_new update_position_and_type!(dt.defectsN[assignment_antivortices[i]],locN_new[i],typeN_new[i],zoom(thetas,lattice,locN_new[i]...,WINDOW)[2]) end

    elseif N_new > 0 && N_old == 0
        for i in 1:Np_new add_defect!(dt,charge=chargeP_new[i],type=typeP_new[i],loc=locP_new[i],thetas_zoom=zoom(thetas,lattice,locP_new[i]...,WINDOW)[2]) end
        for i in 1:Nn_new add_defect!(dt,charge=chargeN_new[i],type=typeN_new[i],loc=locN_new[i],thetas_zoom=zoom(thetas,lattice,locN_new[i]...,WINDOW)[2]) end

    elseif N_new == 0 && N_old > 0 # (+)(-) >> plus rien
        id_just_annihilated_defectP,id_just_annihilated_defectM = ID_active_defects(dt) # seek for not yet annihilated defects

        for i in id_just_annihilated_defectP dt.defectsP[i].annihilation_time = dt.current_time end
        for i in id_just_annihilated_defectM dt.defectsN[i].annihilation_time = dt.current_time end
        # annihilation_time is already taken care of in the annihilate_defects function
        dt = annihilate_defects(dt::DefectTracker,id_just_annihilated_defectP,lattice)

    elseif Np_new > 0 && Np_old > 0 && Nn_old > 0 && Nn_new == 0  # (+)(+)(-) >> (+) par exemple
        assignment_vortices = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_vortices) update_position_and_type!(dt.defectsP[assignment_vortices[i]],locP_new[i],typeP_new[i],zoom(thetas,lattice,locP_new[i]...,WINDOW)[2]) end
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
        dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)

    elseif Nn_new > 0 && Nn_old > 0 && Np_old > 0 && Np_new == 0  # (+)(-)(-) >> (-) par exemple
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")
        # Update living antivortices. NB : the annihilated antivortex is absent from the assignment vector : proceed without the condition "≠ 0"
        for i in eachindex(assignment_antivortices) update_position_and_type!(dt.defectsN[assignment_antivortices[i]],locN_new[i],typeN_new[i],zoom(thetas,lattice,locN_new[i]...,WINDOW)[2]) end
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

        dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)
    else # end of special cases

    # GENERAL TREATMENT
        assignment_vortices     = pair_up_hungarian(dt,locP_new,locP_old,lattice,"+")
        assignment_antivortices = pair_up_hungarian(dt,locN_new,locN_old,lattice,"-")

        # CASE 1 : no creation, no annihilation : simply update the data structure
        if N_new == N_old
            for i in 1:Np_new update_position_and_type!(dt.defectsP[assignment_vortices[i]],locP_new[i],typeP_new[i],zoom(thetas,lattice,locP_new[i]...,WINDOW)[2]) end
            for i in 1:Nn_new update_position_and_type!(dt.defectsN[assignment_antivortices[i]],locN_new[i],typeN_new[i],zoom(thetas,lattice,locN_new[i]...,WINDOW)[2]) end

        # CASE 2 : creation !
    elseif N_new > N_old
            # Take care of the newly created defects
            ind_created_vortex = findall(iszero,assignment_vortices) # newly created vortex -> the assignment vector contains a 0
            loc_created_vortex = vortices_new[ind_created_vortex]
            for j in each(loc_created_vortex) add_defect!(dt,charge=chargeP_new[j],type=typeP_new[j],loc=loc_created_vortex[j][1:2],thetas_zoom=zoom(thetas,lattice,loc_created_vortex[j][1:2]...,WINDOW)[2]) end

            ind_created_antivortex = findall(iszero,assignment_antivortices)
            loc_created_antivortex = antivortices_new[ind_created_antivortex]
            for j in each(loc_created_antivortex) add_defect!(dt,charge=chargeN_new[j],type=typeN_new[j],loc=loc_created_antivortex[j][1:2],thetas_zoom=zoom(thetas,lattice,loc_created_antivortex[j][1:2]...,WINDOW)[2]) end

            # Update the ancient defects' positions
            for i in eachindex(assignment_vortices)
                if assignment_vortices[i] ≠ 0 # avoid newly created defects
                    update_position_and_type!(dt.defectsP[assignment_vortices[i]],locP_new[i],typeP_new[i],zoom(thetas,lattice,locP_new[i]...,WINDOW)[2])
                end
            end
            for i in eachindex(assignment_antivortices)
                if assignment_antivortices[i] ≠ 0 # avoid newly created defects
                    update_position_and_type!(dt.defectsN[assignment_antivortices[i]],locN_new[i],typeN_new[i],zoom(thetas,lattice,locN_new[i]...,WINDOW)[2])
                end
            end

        # CASE 3 : annihilation !
    elseif N_new < N_old
             # Update living vortices. NB : the annihilated vortex is absent from the assignment vector : proceed without the condition "≠ 0"
             for i in eachindex(assignment_vortices)     update_position_and_type!(dt.defectsP[assignment_vortices[i]],locP_new[i],typeP_new[i],zoom(thetas,lattice,locP_new[i]...,WINDOW)[2]) end
             for i in eachindex(assignment_antivortices) update_position_and_type!(dt.defectsN[assignment_antivortices[i]],locN_new[i],typeN_new[i],zoom(thetas,lattice,locN_new[i]...,WINDOW)[2]) end

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
            dt = annihilate_defects(dt,ID_annihilated_vortices,lattice)
        end # end of general treatment
    end # end of special cases & general treatment
    return dt
end

function MSD(dfts::Vector{Union{Missing,DefectTracker}},lattice::Abstract2DLattice)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    maxlength = maximum([maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)]) for dft in dfts[indices]])
    MSD_P   = NaN*zeros(length(indices),maxlength)
    MSD_N   = NaN*zeros(length(indices),maxlength)
    MSD_all = NaN*zeros(length(indices),maxlength)
    for i in 1:length(indices)
        msd_all, msd_p, msd_n = MSD(dfts[indices[i]],lattice)
        MSD_P[i,1:length(msd_p)] = msd_p
        MSD_N[i,1:length(msd_n)] = msd_n
        MSD_all[i,1:length(msd_all)] = msd_all
    end

    MSD_P_avg = nanmean(MSD_P,1)[1,:]
    MSD_N_avg = nanmean(MSD_N,1)[1,:]
    MSD_all_avg = nanmean(MSD_all,1)[1,:]

    return MSD_all_avg,MSD_P_avg,MSD_N_avg
end

function MSD(dft::DefectTracker,lattice::Abstract2DLattice,maxlength=nothing)
    nP = number_defectsP(dft)
    nN = number_defectsN(dft)
    # tmin,tmax = t_bounds(dft) # (tmin,tmax) = timestamps of (first defect creation , last defect annihilation)

    # hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1
    if isnothing(maxlength)
        maxlength = maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)])
    end
    # Compute the SD
    SD_P = NaN*zeros(nP,maxlength)
    SD_N = NaN*zeros(nN,maxlength)
    for n in 1:nP
        defect = dft.defectsP[n]
        tmp = square_displacement(defect,lattice)
        SD_P[n,1:length(tmp)] = tmp
    end
    for n in 1:nN
        defect = dft.defectsN[n]
        tmp = square_displacement(defect,lattice)
        SD_N[n,1:length(tmp)] = tmp
    end

    # Now average to compute the MSD
    MSD_P = nanmean(SD_P,1)[1,:]
    MSD_N = nanmean(SD_N,1)[1,:]
    MSD_all = nanmean(hcat(MSD_P,MSD_N),2)[:]

    return MSD_all, MSD_P, MSD_N
end

function square_displacement(d::Defect,lattice::Abstract2DLattice)
    loc_t0 = first_loc(d)
    return [dist(lattice,loc,loc_t0) for loc in d.pos] .^ 2
end

function interdefect_distance(dft::DefectTracker,defect1::Defect,defect2::Defect,lattice::Abstract2DLattice)
    # TODO take care of case with creation and/or annihilation time different.
    # So far, this care is left to the user...
    # @assert defect1.creation_time == defect2.creation_time
    # @assert defect1.annihilation_time == defect2.annihilation_time
    tmax = min(length(defect1.pos),length(defect2.pos))
    R = [dist(lattice,defect1.pos[t],defect2.pos[t]) for t in 1:tmax]
    return R
end

function mean_distance_to_annihilator(dfts::Vector{Union{Missing,DefectTracker}},lattice::Abstract2DLattice)
    indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
    for i in 1:length(dfts)
        if !ismissing(dfts[i]) push!(indices,i) end
    end
    Rs = [mean_distance_to_annihilator(dfts[indices[n]],lattice) for n in 1:length(indices)]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1]
end

function mean_distance_to_annihilator(dft::DefectTracker,lattice::Abstract2DLattice)
    nP = number_defectsP(dft)
    Rs = [distance_to_annihilator(dft,n,lattice) for n in 1:nP]
    Rs_matrix = vector_of_vector2matrix(Rs)
    return nanmean(Rs_matrix,2)[:,1]
end

function distance_to_annihilator(dft::DefectTracker,id1::Int,lattice::Abstract2DLattice;reversed=true)
    R = interdefect_distance(dft,dft.defectsP[id1],dft.defectsN[dft.defectsP[id1].id_annihilator],lattice)
    if reversed reverse!(R) end
    return R
end

## Small helpful methods for scripts
function number_defects(model::AbstractModel,lattice::Abstract2DLattice)
    a,b = spot_defects(thetas,T,BC)
    return length(a) + length(b)
end

function theta_mid(x::T,y::T,symm)::T  where T<:AbstractFloat
    arcl = arclength(x,y,symm)
    return x + arcl/2 , y - arcl/2
end

## Infer µ
function infer_mu(thetas;q,window=WINDOW,decay=true)
    if decay && q > 0 return infer_mu_decay(thetas,q=q,window=window)
    else return infer_mu_0(thetas,q=q,window=window)
    end
end
function infer_mu_0(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    # tmp = Complex(0)
    range = 2:L-1
    for j in range, i in range
        muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
        # i<->j irrelevant because i,j and j,i have the same weight for "mean" operation
        # tmp += exp(im*muss[i,j] - 0.25sqrt((i-window)^2 + (j-window)^2))
    end
    # muss[window,window] = 0
    corrmat = zeros(size(thetas)) ; corrmat[window,window] = 1
    moyenne = angle(mean(exp.(im*muss[range,range])-0corrmat[range,range]))
    if     abs(q) == 1   correction = pi - 0.33228605 # 0.2 works perfectly with corrmat, 0.33228605 was the original cst
    elseif abs(q) == 1/2 correction = 0.8
    end
    #= To be honest, I don't know where the shifts might come from,
    In the beggining, I thought maybe from the spin at the center of the defect,
    where theta should not be defined. But if one changes mean->sum and adds the condition
    "if (i == j == window)", the result only becomes weirder... =#
    # return mod.(muss,2pi)
    return mod(moyenne .+ correction,2π)
end
function infer_mu_decay(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    tmp = Complex(0)
    range = 2:L-1
    for j in range, i in range
        muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
        # i<->j irrelevant because i,j and j,i have the same weight for "mean" or "sum" operation
        tmp += exp(im*muss[i,j] - sqrt((i-window)^2 + (j-window)^2))
    end
    moyenne = angle(tmp)
    if     q == 1   correction = pi - 0.603228605 # 0.33228605 was the original cst
    elseif q == 1/2 correction = 0.8
    elseif q == -1/2 correction = 0.25
    end
    #= To be honest, I don't know where the shifts might come from,
    In the beggining, I thought maybe from the spin at the center of the defect,
    where theta should not be defined. But if one changes mean->sum and adds the condition
    "if !(i == j == window)", the result only becomes weirder... =#
    return mod(moyenne .+ correction,2π)
end
