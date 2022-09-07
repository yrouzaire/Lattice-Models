include("lattices.jl");
include("models.jl");
include("core_methods.jl");
import StatsBase.sample
import Plots.@animate

## ------------------------ Initializations  ------------------------
function init_thetas(model::AbstractModel,space::AbstractLattice,params)
    @unpack L,init,q,r0,type1defect,type2defect = params
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
        thetas = create_pair_vortices(L,r=r,q=q,type=type2defect)
    elseif init in ["2pairs" , "2pair"]
        thetas = create_2pairs_vortices(L,r=r,q=q,type=type2defect)
    else error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\",\"isolated\" , \"pair\" , \"2pair\" or \"lowtemp_nematic/nematic_order\" .")
    end
    if isa(model,MovingXY) make_holes!(thetas,model.rho) end
    return float_type.(thetas)
end

function create_single_defect(L,x0=round(Int,L/2),y0=round(Int,L/2);q=1,type="source")
    @assert (q > 0 && type in ["source","sink","clockwise","counterclockwise"]) || (q < 0 && type in ["convergent","divergent","threefold1","threefold2"])
    thetas = zeros(L,L)
    for y in 1:L , x in 1:L
        # q > 0
        if     type == "clockwise"            thetas[x,y] = q * atan(y-y0,x-x0)
        elseif type == "counterclockwise"     thetas[x,y] = q * atan(y-y0,x-x0) + pi
        elseif type == "sink"                 thetas[x,y] = q * atan(y-y0,x-x0) + pi/2
        elseif type == "source"               thetas[x,y] = q * atan(y-y0,x-x0) - pi/2
        # q = -1/2 (when q = -1, I think they are all equivalent)
        elseif type == "threefold1"           thetas[x,y] = q * atan(y-y0,x-x0)
        elseif type == "threefold2"           thetas[x,y] = q * atan(y-y0,x-x0) + pi
        elseif type == "divergent"            thetas[x,y] = q * atan(y-y0,x-x0) + pi/2
        elseif type == "convergent"           thetas[x,y] = q * atan(y-y0,x-x0) - pi/2
        end
    end
    return thetas
end

function create_pair_vortices(L;r0=Int(L/2),q,type)
    #= Check for meaningfulness of the defaults separation,
    otherwise the defaults will annihilate with relaxation =#
    @assert r0 ≤ 0.5L  "Error : r0 > L/2. "
    @assert iseven(r0) "Error : r0 has to be even. "

    thetas = create_single_defect(L,round(Int,L/2+r0/2),round(Int,L/2),q=+q,type=type[1]) +
             create_single_defect(L,round(Int,L/2-r0/2),round(Int,L/2),q=-q,type=type[2])
    return thetas
end

function create_2pairs_vortices(L;r0,q,type)
    @assert isinteger(r0/4) "R0 has to be a multiple of 4 !"
    return create_pair_vortices(L;r0,q,type) + create_pair_vortices(L;r0=r0/2,q,type)'
end

function make_holes!(thetas,rho)
    N = length(thetas)
    holes = sample(1:N,round(Int,(1-rho)*N),replace=false)
    thetas[holes] .= NaN
    return thetas
end

## ------------------------ Visualization  ------------------------
function plot_theta(thetas::Matrix{<:AbstractFloat},model::AbstractModel,lattice::AbstractLattice;defects=false,title="",colorbar=true,cols = cgrad([:black,:blue,:green,:orange,:red,:black]),size=(485,400))
    symm = sym(model)
    p = heatmap(mod.(thetas',symm),c=cols,clims=(0,symm),size=size,
        colorbar=colorbar,colorbartitle="θ",title=title,aspect_ratio=1)

    if defects
        defects_p,defects_m = spot_defects(thetas,model,lattice)
        highlight_defects!(p,lattice.L,defects_p,defects_m)
    end
    xlims!((0,lattice.L))
    ylims!((0,lattice.L))
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

function display_quiver!(p,thetas,window)
    p
    for j in 1:2window+1
        quiver!(j*ones(2window+1),collect(1:2window+1),
        quiver=(cos.(thetas'[j,:]),-sin.(thetas'[j,:])),
        c=:white,lw=0.8)
    end
    return p
end
# Note : after returning the plots with quiver, one has to add xlims!(1,2window+1) ; ylims!(1,2window+1)

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
