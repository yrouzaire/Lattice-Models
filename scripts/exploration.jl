using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Tracking a pair over time
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = SquareLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
plot_thetas(thetas,model,lattice,defects=true)
dft = DefectTracker(thetas,model,lattice)
## Simple simulation and plots from hightemp
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
times = [200,300,500,800,1000]
times = [5,10,20,30,50,80,100]
# @elapsed z = update!(thetas,model,lattice,100)
token = 1
while model.t < times[end]
    update!(thetas,model,lattice,100)
    if model.t ≥ times[token]
        token = min(token+1,length(times))
        display(plot_thetas(thetas,model,lattice,defects=true))
    end
end
@elapsed z = update!(thetas,model,lattice,410)
plot_thetas(thetas,model,lattice)
spot_defects(thetas,model,lattice)
plot_thetas(thetas,model,lattice,defects=true)
