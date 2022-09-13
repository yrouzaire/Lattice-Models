using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

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

## Parameters
include(srcdir("../parameters.jl"));

params["algo"] = "A"
    params["rho"] = 0.95
    params["width_proposal"] = 0.1
    params["symmetry"] = "nematic"
    params["A"] = 3

    params["init"] = "hightemp"
    params["type1defect"] = "source"
    params["type2defect"] = "pair4"
    params["q"] = 1/2
    params["r0"] = 40

model = XY(params)
    lattice = SquareLattice(L,periodic=true)
    thetas = init_thetas(lattice,params=params)
update!(thetas,model,lattice,1000)
    plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,2000)
    plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,3000)
    plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,4000)
    plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,5000)
    plot_thetas(thetas,model,lattice,defects=true)

spot_defects(thetas,model,lattice)


## Do some movies
lattice = TriangularLattice(L)
saving_times = vcat(0:10:700,700:100:20000) ; transients = saving_times[end]/3

# params["rho"] = 1
#     model = MCXY(params)
#     thetas = init_thetas(lattice,params=params)
#     z = @elapsed anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
#     mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MCXY_rho1.mp4")
#     prinz(z)
#
# params["rho"] = 0.95
#     model = MCXY(params)
#     thetas = init_thetas(lattice,params=params)
#     anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
#     mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MCXY_rho0.95.mp4")

params["rho"] = 1
    params["A"] = 0.5
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    z = @elapsed anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")
    prinz(z)
    
    params["A"] = 0
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")


    params["A"] = 1
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")

params["rho"] = 0.95
    params["A"] = 0
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")

    params["A"] = 0.5
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")

    params["A"] = 1
    model = MCXY(params)
    thetas = init_thetas(lattice,params=params)
    anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
    mp4(anim,"D:/Documents/Research/projects/LatticeModels/films/MovXY_rho$(params["rho"])_A$(params["A"]).mp4")
