cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

## SPP
params["symmetry"] = "nematic" ; params["rho"] = 0.95 ; params["A"] = 2
model = SPP(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice,tmax=500)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)

tmax = 1000
saving_times = vcat(model.t:1:tmax) ; transients = Inf #saving_times[end]/3
anim = @animate for t in saving_times
    println("$(round(t/tmax*100,digits=2)) %")
    update!(thetas,model,lattice,tmax=t)  # updates until time = t
    if t<transients #p = plot_thetas(thetas,model,lattice,defects=false,size=(512,512))
        zoom_quiver(thetas,model,lattice,93,120,12)
    else            p = plot_thetas(thetas,model,lattice,defects=defects,size=(512,512))

    end
end
film = mp4(anim,"films/movies_and_analysis/test2_zoom_ralenti.mp4",fps=5)
p = plot_thetas(thetas,model,lattice,defects=false)
zoom_quiver(thetas,model,lattice,7,145)

update!(thetas,model,lattice)
    zoom_quiver(thetas,model,lattice,92,120,12)
