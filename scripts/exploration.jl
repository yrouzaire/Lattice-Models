cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Movies
include(srcdir("../parameters.jl"));

model = MovingXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
saving_times = 0:100:20000 ; transients = saving_times[end]/2
z = @elapsed anim_extensile = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
prinz(z)
mp4(anim_extensile,datadir("../films/active_extensile_A2.mp4"))


## Verifions que les défauts soient bien les bons
global const WINDOW = 7
using Flux:onecold, Chain, Dense, softmax
global const NN_positive = load("NN_positive_12_defects_N1000_W7.jld2","NN")
global const possible_positive_defects = load("NN_positive_12_defects_N1000_W7.jld2","possible_defects")
global const NN_negative = load("NN_negative_12_defects_N1000_W7.jld2","NN")
global const possible_negative_defects = load("NN_negative_12_defects_N1000_W7.jld2","possible_defects")

include(srcdir("../parameters.jl"));
    model = MovingXY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)

plot_thetas(thetas,model,lattice)
update!(thetas,model,lattice,Int(2E3))
    plot_thetas(thetas,model,lattice)
plot_thetas(thetas,model,lattice,defects=true)

dft = DefectTracker(thetas,model,lattice,find_type=true)
number_defects_types(dft)

n = 1
n += 1
    d = dft.defectsP[n]
    zoom_quiver(thetas,model,lattice,last_loc(d)...,11)
    thetas_zoom = zoom(thetas,lattice,last_loc(d)...,WINDOW)[2]
    preconditionning!(thetas_zoom,model,lattice)
    pred = NN_positive(vec(thetas_zoom))
    confiance = maximum(pred)/sum(pred)
    title!(possible_positive_defects[onecold(pred)]*" confiance: $(round(confiance,digits=2))")

# function get_types(dft::DefectTracker)
#
# end


## Tracking defects over time (first single, then pair, then hightemp)
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
update!(thetas,model,lattice,20)
p=plot_thetas((thetas),model,lattice,defects=false)
    # display_quiver!(p,(thetas),13)

dft = DefectTracker(thetas,model,lattice)
update_DefectTracker!(dft,thetas,model,lattice)
number_defects_types(dft)
# dft.defectsP[1].type[1]
# dft.defectsN[1].type[1]

update_and_track!(thetas,model,lattice,dft,4000,50)
dft.defectsP[7].type
dft.defectsN[1].type
number_defects_types(dft)
[last_loc(dt.defectsP[i])]


zoom_quiver(thetas,model,lattice,last_loc(dft.defectsP[1])...,7)
zoom_quiver(thetas,model,lattice,last_loc(dft.defectsN[1])...,7)
zoomN = zoom(thetas,lattice,last_loc(dft.defectsN[1])...,7)[2]
ind_type = onecold(NN(vec(zoomN)))
possible_defects[ind_type]




tmax,every = 100,10
z = @elapsed update_and_track_plot!(thetas,model,lattice,dft,tmax,every)
prinz(z)
dft

## Simple simulation and plots from hightemp
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(model,lattice,params_init=params_init)
    times = [5,10,20,30,50,80,100,150,200,250,300,350,400,450,500]
    # times = [5,10,20,30]
    # @elapsed z = update!(thetas,model,lattice,100)
    token = 1
    while model.t < times[end]
        update!(thetas,model,lattice)
        if model.t ≥ times[token]
            token = min(token+1,length(times))
            # display(plot_thetas(thetas,model,lattice,defects=true))
            display(zoom_quiver(thetas,model,lattice,spot_defects(thetas,model,lattice)[1][2][1:2]...,9))
        end
    end

@elapsed z = update!(thetas,model,lattice,10)
plot_thetas(thetas,model,lattice)
spot_defects(thetas,model,lattice)
plot_thetas(thetas,model,lattice,defects=true)
