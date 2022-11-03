using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

#= First Test: create an isolated defect and track it over time
to see what happen to µ(t) in controlled conditions =#

lattice = TriangularLattice(L)
params["symmetry"] = "nematic" ; params["rho"] = 1 ; params["A"] = 0
    params_init["init"] = "single" ; params_init["type1defect"] = 1 # initial µ
    # params_init["init"] = "pair" ; params_init["type2defect"] = [0,pi/2] # initial µ
    model = LangevinXY(params)
    thetas = init_thetas(model,lattice,params_init=params_init)
    # p = plot_thetas(thetas,model,lattice,defects=false)

    dft = DefectTracker(thetas,model,lattice,find_type=true)
    update_and_track!(thetas,model,lattice,dft,25,.5,find_type=true)
    plot!(dft.defectsP[1].type,m=false)
dft.defectsP[1].type
dft.defectsP[1].pos
zoom_quiver(thetas,model,lattice,65,70,9)

update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,1000)  # updates until time = t
model.t


histogram(last_types(dft),bins=50,normalize=true)
dft

~,thetas_zoom= zoom(thetas,lattice,64.5,68.8)
zoom_quiver(thetas_zoom,model,lattice,8,8)
recon = reshape(DAE_positive12(reshape(thetas_zoom,(15,15,1,1))),(15,15))
zoom_quiver(recon,model,lattice,8,8)
infer_mu(recon,1/2)
