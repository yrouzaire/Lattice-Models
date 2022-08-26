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
