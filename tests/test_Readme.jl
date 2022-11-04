using DrWatson ; @quickactivate "LatticeModels" # tell Julia to load the correct package versions
include(srcdir("LatticeModels.jl")) # loads the entire code, might precompile some packages

using Plots,ColorSchemes,LaTeXStrings
# Load the plotting package. Can be very slow first time it is run, could take up to 30s - 1min
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

#  Declare your lattice: SquareLattice(L) or TriangularLattice(L).
lattice = SquareLattice(L)

# Declare your model: LangevinXY(params), MonteCarloXY(params), VisionXY(params), SPP(params) etc
model = LangevinXY(params)

# Initialisation of the theta field:
thetas = init_thetas(model,lattice,params_init=params_init)

model.t # Now the current time of the model = 0
# Update the model:
update!(thetas,model,lattice) # For one time step only
model.t # Now the current time of the model = dt
duration = 5
update!(thetas,model,lattice,duration) # For a given duration, here 1:
model.t # Now the current time of the model = dt + duration
tmaxx = 50
update!(thetas,model,lattice,tmax=tmaxx) # Until a given time, here 50
model.t # Now the current time of the model = tmaxx, no matter what was done before


# Visualise the system
plot_thetas(thetas,model,lattice,defects=false)
plot_thetas(thetas,model,lattice,defects=true) # circle = +1 defect, triangle = -1 defect

# Some measurement on the system
polarOP, nematicOp = OP(thetas) # polar and nematic order parameters
C = corr(thetas,model,lattice) # the spatial correlation function C(r)
threshold = exp(-1)
ξ = corr_length(C,seuil=threshold) # the correlation length ξ (xi)
plot(xlabel="r",ylabel="C(r)")
    plot!(1:L/2,C,axis=:log)
    hline!([threshold],c=:grey,lw=0.5)
    vline!([ξ],c=:grey,lw=0.5,line=:dash)
    annotate!((1.5,1.1threshold,text("y = threshold",7,:left)))
    annotate!((8,0.1,text("r = ξ",7,:center,90.)))

# Where are the defects ?
positively_charged_defects,negatively_charged_defects = spot_defects(thetas,model,lattice)
positively_charged_defects[1] # charge, x, y, type (don't care about it)
