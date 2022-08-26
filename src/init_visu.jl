include("lattice_general.jl");
include("models.jl");
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
plot()

## ------------------------ Tests  ------------------------
L = 200
model1 = XY(L,0.1,"polar")
lattice = TriangularLattice(L)
init_thetas!(model1,lattice,"hightemp")

## ------------------------ Initializations  ------------------------
function init_thetas!(model::AbstractModel,space::AbstractLattice,init::String)
    @assert model.L == space.L
    L = model.L
    # model.symmetry == "polar" ? symm = 2 : symm = 1
    if init == "hightemp" || "disorder"
        model.thetas = 2π*rand(L,L)
    elseif init == "lowtemp" || "polar_order"
        model.thetas = zeros(L,L)
    elseif init == "lowtemp_nematic" || "nematic_order"
        model.thetas = rand(Bool,L,L)*π
    # elseif init == "isolated"  thetas = create_single_defect(args...,L,round(Int,L/2),round(Int,L/2)) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
    # elseif init == "pair"      thetas = create_pair_vortices(L,r0,args...)
    # elseif init == "2pair"     thetas = create_pair_vortices(L,r0,args...) + create_pair_vortices(L,round(Int,r0/2),args...)'
    else error("ERROR : Type of initialisation unknown. Choose among \"hightemp/order\",\"lowtemp/polar_order\",\"isolated\" , \"pair\" , \"2pair\" or \"lowtemp_nematic/nematic_order\" .")
    end
    return nothing
end

## ------------------------ Visualization  ------------------------
function plot_theta(model::AbstractModel,space::AbstractLattice;defects=false,title="",colorbar=true,cols = cgrad([:black,:blue,:green,:orange,:red,:black]))
    model.symmetry == "polar" ? symm = 2π : symm = π
    p = heatmap(mod.(model.thetas',symm),c=cols,clims=(0,symm),size=(512,512),
        colorbar=colorbar,title=title,aspect_ratio=1)

    if defects
        defects_p,defects_m = spot_defects(model.thetas,model.T,space.periodic)
        highlight_defects!(p,model.L,defects_p,defects_m)
    end

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
