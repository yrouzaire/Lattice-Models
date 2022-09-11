using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Parameters
include(srcdir("../parameters.jl"));

#= Qualitative tests that the algo has to pass
1. Numerical Stability from lowtemp
2. Numerical Stability from lowtemp_nematic
3. Good Looking for XY limit from hightemp
4. Good Looking for nematic defects (comet shaped and triangle)
=#

#= Quantitative tests that the algo has to pass
1. XY T_kt=0.45 and C(r,t)
2.
=#

#= Summary of Algorithms
A. Upon collision, align nematically wrt all NN.
B. Upon collision, align nematically wrt collided spin
C. Upon collision, align F/AF wrt collided spin
In all cases,
    - only the colliding spin is updated.
    - the proposal is taken from a gaussian centered on the current orientation
    - the proposal is accepted with Metropolis proba (to ensure a well controlled limit at eq)
=#

params["algo"] = "A"
params["width_proposal"] = 0.05
params["symmetry"] = "nematic"

params["init"] = "hightemp"
    params["A"] = 1.
    params["q"] = -1/2

model = MovingXY(params)
    lattice = TriangularLattice(L,periodic=true)
    thetas = init_thetas(lattice,params=params)
update!(thetas,model,lattice,3000)
    plot_thetas(thetas,model,lattice,defects=true)
