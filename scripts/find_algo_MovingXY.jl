using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Parameters
include(srcdir("../parameters.jl"));

#= Tests that the algo has to pass
1. Numerical Stability from lowtemp
2. Numerical Stability from lowtemp_nematic
3. Good Looking for XY limit from hightemp
4. Good Looking for nematic defects (comet shaped and triangle)
=#

## Test 1 Numerical Stability from lowtemp
params["init"] = "lowtemp"
    params["antiferro"] = false
    params["q"] = -1/2
    params["width_proposal"] = 0.1
    params["symmetry"] = "nematic"

model = XY(params)
    lattice = TriangularLattice(L,periodic=true)
    thetas = init_thetas(lattice,params=params)
    update!(thetas,model,lattice,1)
    plot_thetas(thetas,model,lattice,defects=true)
