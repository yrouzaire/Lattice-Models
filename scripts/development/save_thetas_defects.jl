using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

target_number_defects = [20,14,10,6,4,2]
configurations = Vector{Matrix{<:AbstractFloat}}(undef,length(target_number_defects))
every = 1.0

## XY Model
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    nb_def = Inf
    while nb_def > target_number_defects[end]
        update!(thetas,model,lattice,model.t + every)
        nb_def = number_defects(thetas,model,lattice)
