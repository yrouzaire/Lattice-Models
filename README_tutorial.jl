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
#= WARNING: plotting the defects can take enormous time if they are too many.
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
every = 1 ; tmax = 20 ; transients = 50 # defects are not plotted before t ≥ transients (if defects=true)
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

params_init["init"] = "lowtemp" # equivalent of "ordered"
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)

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

## Run different simulations and save them
using JLD2 # saving package, will save data under .jld2 format
include(srcdir("../parameters.jl"));
R = 1 # number of indep realisations
tmax = 1000
# every  = 2 ; saving_times = every:every:tmax # lin spaced
nb_sav_pts  = 8 ; saving_times = logspace(1,tmax,nb_sav_pts) # log spaced
magnetisations = zeros(length(saving_times),R)
corr_functions = zeros(Int(L/2),length(saving_times),R)
corr_lengths   = zeros(length(saving_times),R)

z = @elapsed for r in 1:R
    println("Simulation $r/$R started.")
    # Prepare your setup
    model = LangevinXY(params) # sets t = 0
    lattice = SquareLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    token = 1
    while model.t < tmax # run the actual simulation
        update!(thetas,model,lattice)
        if model.t ≥ saving_times[token] # once you hit a saving time, save quantites of interest
            println("t=$(round(model.t,digits=1)), $(round(model.t/tmax*100,digits=1))% achieved.")
            magnetisations[token,r]   = OP(thetas)[1] # polar Order Parameter
            corr_functions[:,token,r] = corr(thetas,model,lattice) # correlation function C(r)
            corr_lengths[token,r]     = corr_length(corr_functions[:,token,r]) # correlation length computed from C(r)
            token = min(token+1,length(saving_times)) # next saving_time
        end
    end
end
prinz(z) # takes 4.5 mn on my machine (i7,16Go RAM) for L = 100, R = 1, tmax = 1000
filename = datadir("myfirstdata.jld2")
JLD2.@save filename magnetisations corr_functions corr_lengths params # save data
JLD2.@load filename magnetisations corr_functions corr_lengths params # load data
corr_functions_avg = mean(corr_functions,dims=3)
magnetisations_avg = mean(magnetisations,dims=2)
corrlength_avg = mean(corr_lengths,dims=2)

p=plot(xlabel="r",ylabel="C(r,t)",axis=:log,ylims=(1E-2,1.2))
    for tt in 1:length(saving_times)
        plot!(1:params["L"]/2,remove_negative(corr_functions_avg[:,tt,1]))
    end
    p # displays the figure "p"

plot(saving_times,magnetisations_avg,axis=:log,m=true,xlabel="t",ylabel="Magnetisation P(t)")
plot(saving_times,corrlength_avg,axis=:log,m=true,xlabel="t",ylabel="ξ(t)")



corr_length(corr_functions_avg[:,end,1])
