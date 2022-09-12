include("lattices.jl");
include("models.jl");
include("core_methods.jl");
import StatsBase.sample
import Plots.@animate

## ------------------------ Initializations  ------------------------
function init_thetas(space;params)
    @unpack L,rho,init,q,r0,type1defect,type2defect,float_type = params
    if init in ["hightemp" , "disorder"]
        thetas = 2π*rand(L,L)
    elseif init in ["lowtemp" , "polar_order"]
        thetas = zeros(L,L)
    elseif init in ["lowtemp_nematic" , "nematic_order"]
        thetas = rand(Bool,L,L)*π
    elseif init in ["isolated" , "single"]
        thetas = create_single_defect(L,round(Int,L/2),round(Int,L/2),q=q,type=type1defect) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
        space.periodic = false
    elseif init == "pair"
        thetas = create_pair_vortices(L,r0=r0,q=q,type=type2defect)
    elseif init in ["2pairs" , "2pair"]
        thetas = create_2pairs_vortices(L,r0=r0,q=q,type=type2defect)
    else error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\",\"isolated\" , \"pair\" , \"2pair\" or \"lowtemp_nematic/nematic_order\" .")
    end
    if model.rho < 1 make_holes!(thetas,rho) end
    return float_type.(thetas)
end

function create_single_defect(L,x0=round(Int,L/2),y0=round(Int,L/2);q=1,type)
    condition1 = (q > 0 && type in ["source","sink","clockwise","counterclockwise"])
    condition2 = (q < 0 && type in ["split","join","convergent","divergent","threefold1","threefold2"])
    @assert condition1 || condition2

    thetas = zeros(L,L)
    for y in 1:L , x in 1:L
        # q > 0
        if     type == "counterclockwise"     thetas[x,y] = q * atan(y-y0,x-x0)
        elseif type == "clockwise"            thetas[x,y] = q * atan(y-y0,x-x0) + pi
        elseif type == "sink"                 thetas[x,y] = q * atan(y-y0,x-x0) + pi/2
        elseif type == "source"               thetas[x,y] = q * atan(y-y0,x-x0) - pi/2
        # q = -1/2 (when q = -1, I think they are all equivalent)
        elseif type in ["threefold2"]          thetas[x,y] = q * atan(y-y0,x-x0)
        elseif type in ["threefold1"]          thetas[x,y] = q * atan(y-y0,x-x0) + pi
        elseif type in ["divergent" ,"split"]     thetas[x,y] = q * atan(y-y0,x-x0) + pi/2
        elseif type in ["convergent","join"]      thetas[x,y] = q * atan(y-y0,x-x0) - pi/2
        end
    end
    return thetas
end

function create_pair_vortices(L;r0=Int(L/2),q,type)
    #= Check for meaningfulness of the defaults separation,
    otherwise the defaults will annihilate with relaxation =#
    @assert r0 ≤ 0.5L  "Error : r0 > L/2. "
    @assert iseven(r0) "Error : r0 has to be even. "

    if isa(type,String)
        if     type in ["pair1","source_split"] type_pos,type_neg = "source","threefold1"
        elseif type in ["pair2","sink_join"]    type_pos,type_neg = "source","threefold2"
        elseif type in ["pair3","cw_31"]        type_pos,type_neg = "source","join"
        elseif type in ["pair4","ccw_32"]       type_pos,type_neg = "source","split"
        else error("Type Unknown!")
        end
    elseif isa(type,Vector{String})
        type_pos , type_neg = type
    else error("Type Unknown!")
    end

    thetas = create_single_defect(L,round(Int,L/2+r0/2),round(Int,L/2),q=+q,type=type_pos) +
             create_single_defect(L,round(Int,L/2-r0/2),round(Int,L/2),q=-q,type=type_neg)

    # return smooth_border!(thetas)
    return (thetas)
end

function smooth_border!(thetas)
    cst = round(Int,0.02*size(thetas,1)) # has to be even

    for i in 1:L
        sample1 = thetas[i,cst]
        sample2 = thetas[i,L-cst]
        arcl = arclength(sample1,sample2,2pi) # Oddly, pi or 2pi is not important even when symmetry is polar or nematic
        smoothed_values = [sample1 + arcl*n/(2cst+1) for n in 1:(2cst+1)]
        for j in 1:2cst+1
            thetas[i,mod1(Int(cst-j+1),L)] = smoothed_values[j]
        end
    end
    return thetas # in fact, this procedure is fairly insensible to relax!() since already smooth
end

function create_2pairs_vortices(L;r0,q,type)
    @assert isinteger(r0/4) "R0 has to be a multiple of 4 !"
    return create_pair_vortices(L,r0=r0,q=q,type=type) + create_pair_vortices(L,r0=r0/2,q=q,type=type)'
end

function make_holes!(thetas,rho)
    N = length(thetas)
    holes = sample(1:N,round(Int,(1-rho)*N),replace=false)
    thetas[holes] .= NaN
    return thetas
end

## ------------------------ Visualization  ------------------------
function plot_thetas(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::AbstractLattice;defects=false,title="",colorbar=true,cols = cgrad([:black,:blue,:green,:orange,:red,:black]),size=(485,400))
    if     model.symmetry == "nematic" modd = pi
    elseif model.symmetry == "polar"   modd = 2pi
    end
    p = heatmap(mod.(thetas',modd),c=cols,clims=(0,modd),size=size,
        colorbar=colorbar,colorbartitle="θ",title=title,aspect_ratio=1)

    if defects
        defects_p,defects_m = spot_defects(thetas,model,lattice)
        locP = [defects_p[i][1:2] for i in each(defects_p)]
        locN = [defects_m[i][1:2] for i in each(defects_m)]
        highlight_defects!(p,lattice.L,locP,locN)
    end
    # xlims!((0,lattice.L))
    # ylims!((0,lattice.L))
    return p
end

function highlight_defects!(p,L,defects_p,defects_m,symbP=:circle,symbM=:utriangle)
    for defect in defects_p
        scatter!((defect), m = (8, 12.0, symbP,:transparent, stroke(1.2, :grey85)))
    end
    for defect in defects_m
        scatter!((defect), m = (8, 12.0, symbM,:transparent, stroke(1.2, :grey85)))
    end
    xlims!((1,L))
    ylims!((1,L))
    return p
end

function display_quiver!(p,thetas_zoom,window)
    p
    for j in 1:2window+1
        quiver!(j*ones(2window+1),collect(1:2window+1),
        # quiver=(sin.(thetas_zoom[j,:]),-cos.(thetas_zoom[j,:])), # original
        quiver=(-sin.(thetas_zoom[j,:]),cos.(thetas_zoom[j,:])),
        c=:white,lw=0.8)
    end
    return p
end
# Note : after returning the plots with quiver, one has to add
# xlims!(1,2window+1) ; ylims!(1,2window+1)

## ------------------------ Movies  ------------------------
function movies(thetas,model,lattice;defects=false,saving_times,transients)
    anim = @animate for t in saving_times
        println("$(round(t/saving_times[end]*100,digits=2)) %")
        update!(thetas,model,lattice,t)  # updates until time = t
        if t<transients p = plot_theta(thetas,model,lattice,defects=false,size=(512,512))
        else            p = plot_theta(thetas,model,lattice,defects=defects,size=(512,512))
        end
    end
    return anim
end
