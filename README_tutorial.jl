using DrWatson ; @quickactivate "LatticeModels" # tells Julia to load the correct package versions
include(srcdir("LatticeModels.jl")) # loads the entire code, might precompile some packages

using Plots,ColorSchemes,LaTeXStrings
# Load the plotting package. Can be quite slow first time it is run, could take up to 30s - 1min
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

#= TODO: Modify parameters.jl with the desired values.
For this tuto, take
    L = 100
    T = 0.1
    symmetry = "polar"
    rho = 1
    dt = 1E-2
    init = "hightemp"
    The rest doesn't matter.
=#
include(srcdir("../parameters.jl")) # Load them

## Declare your lattice: SquareLattice(L) or TriangularLattice(L).
lattice = SquareLattice(L)

## Declare your model: LangevinXY(params), MonteCarloXY(params), VisionXY(params), SPP(params) etc
model = LangevinXY(params)

## Initialisation of the theta field:
thetas = init_thetas(model,lattice,params_init=params_init)

## Temporal evolution
model.t # current time of the model = 0
# Update the model:
update!(thetas,model,lattice) # For one time step only
model.t # Now the current time of the model = dt
#= Indeed, the `update!` function also updates time.
Thus, don't forget to re-instantiate the model if you want to
restart a fresh simulation =#

duration = 15
z = @elapsed update!(thetas,model,lattice,duration) # For a given duration
model.t # Now the current time of the model = dt + duration
prinz(z) # prints the rounded runtime of this function call
tmax = 50
update!(thetas,model,lattice,tmax=tmax) # Until a given time
model.t # Now the current time of the model = tmax, no matter what was done before
# equivalently, if one wants to perform measurements over time
tmax2 = 60
while model.t < tmax2
  update!(thetas,model,lattice)
  # perform measurements
end

# Visualise the system
plot_thetas(thetas,model,lattice,defects=false)
plot_thetas(thetas,model,lattice,defects=true) # circle = +1 defect, triangle = -1 defect
#= WARNING: plotting the defects can take enormous time if they are to many.
It is recommended to first plot with `defects=false` (the default value) to evaluate visually
the number of defects. =#

## Some measurement on the system
polarOP, nematicOp = OP(thetas) # polar and nematic order parameters
C = corr(thetas,model,lattice) # the spatial correlation function C(r)
threshold = exp(-1)
ξ = corr_length(C,seuil=threshold) # the correlation length ξ (xi)
plot(xlabel="r",ylabel="C(r)")
    plot!(1:L/2,remove_negative(C),axis=:log) # remove_negative
    hline!([threshold],c=:grey,lw=0.5)
    vline!([ξ],c=:grey,lw=0.5,line=:dash)
    annotate!((1.5,1.15threshold,text("y = threshold",7,:left)))
    annotate!((0.85ξ,0.1,text("r = ξ",7,:center,90.)))

## Create a movie
lattice = SquareLattice(L)
model = LangevinXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
every = 1 ; tmax = 200 ; transients = 50 # defects are not plotted before t ≥ transients (if defects=true)
saving_times = every:every:tmax
z = @elapsed animation = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
prinz(z) # takes approximately 2 minutes, go grab a coffee
filename = "films/my_first_movie_of_critical_XY_model.mp4"
mp4(animation,filename,fps=20) # creates the file in the given directory
# open manually the .mp4 file

## Visualise the flow field
plot_thetas(thetas,model,lattice)

# At random loc
xlocation, ylocation = (50,50)
half_width_of_window = 10 # not too big because plotting the arrows is damn slow
zoom_quiver(thetas,model,lattice,xlocation,ylocation,half_width_of_window)

# Where are the defects ?
vortices,antivortices = spot_defects(thetas,model,lattice)
nb_defects = length(vortices)

# Around a randomly chosen +1 defect
xlocation, ylocation = vortices[rand(1:nb_defects)][1:2]
    zoom_quiver(thetas,model,lattice,xlocation,ylocation,half_width_of_window)

# Around a randomly chosen -1 defect
xlocation, ylocation = antivortices[rand(1:nb_defects)][1:2]
    zoom_quiver(thetas,model,lattice,xlocation,ylocation,half_width_of_window)

## Explore the different initialisations
params_init["init"] = "hightemp" # equivalent of "disordered"
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
&
params_init["init"] = "lowtemp" # equivalent of "ordered"
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
&
# For a +1 defect, "type1defect" should be "source" or "sink" or "clockwise" or "counterclockwise"
params_init["init"] = "single" ;  params_init["q"] = 1 ; params_init["type1defect"] = "source"
    thetas = init_thetas(model,lattice,params_init=params_init)
    p=plot_thetas(thetas,model,lattice)
zoom_quiver(thetas,model,lattice,50,50)

params_init["init"] = "single" ;  params_init["q"] = +1 ; params_init["type1defect"] = "sink"
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
zoom_quiver(thetas,model,lattice,50,50)

# For a -1 defect, "type1defect" should be "join" or "split" or "threefold1" or "threefold2"
    # For now, it doesn't mak much sense, but it does for some other use of this code.
params_init["init"] = "single" ;  params_init["q"] = -1 ; params_init["type1defect"] = "join"
    thetas = init_thetas(model,lattice,params_init=params_init)
    p=plot_thetas(thetas,model,lattice)
zoom_quiver(thetas,model,lattice,50,50)

params_init["init"] = "single" ;  params_init["q"] = -1 ; params_init["type1defect"] = "threefold1"
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
zoom_quiver(thetas,model,lattice,50,50)

# Pairs of defects
params_init["init"] = "pair" ;  params_init["r0"] = round(Int,L/2) ; params_init["type2defect"] = ["source","join"]
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)

# 2 Pairs
params_init["init"] = "2pair" ;  params_init["r0"] = 40 ; params_init["type2defect"] = ["source","join"]
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
