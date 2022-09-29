cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

#= The goal of this file is to investigate the dynamics and behaviour of
manually created defects (= in vitro) over time, to see which type of defects
actually self propel or not, what is the interaction with holes etc...

Results :
+1 single defect, polar symmetry. All is intuitive.


=#

possible_defects = ["source","sink","clockwise","counterclockwise"]
# possible_defects = ["join","split","threefold1","threefold2"]

## Single defects
include(srcdir("../parameters.jl"));

params["symmetry"] = "nematic"
model = MovingXY(params)
lattice = TriangularLattice(L)
params_init["init"] = "single"
params_init["q"] = 1
params_init["type1defect"] = possible_defects[4]

thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,50)
    plot_thetas(thetas,model,lattice,defects=false)
zoom_quiver(thetas,model,lattice,110,110,15)


## Pairs of defects
include(srcdir("../parameters.jl"));

params["symmetry"] = "nematic"
model = MovingXY(params)
lattice = TriangularLattice(L)
params_init["init"] = "pair"
params_init["q"] = 1/2
params_init["type2defect"] = "pair2"

thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,100)
    plot_thetas(thetas,model,lattice,defects=false)
thetas .-= pi/2 
zoom_quiver(thetas,model,lattice,225,150,15)
zoom_quiver(thetas,model,lattice,75,150,15)
