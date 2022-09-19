cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## With TriangularLattice
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
# plot_thetas(thetas,model,lattice)

# OP(thetas)
update!(thetas,model,lattice,10)
ctri = corr(thetas,model,lattice)
plot(ctri)
corr_length(ctri)

## With SquareLattice
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
# plot_thetas(thetas,model,lattice)

# OP(thetas)
update!(thetas,model,lattice,10)
csqu = corr(thetas,model,lattice)
plot!(csqu)






&
