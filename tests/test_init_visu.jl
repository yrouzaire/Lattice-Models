include("../src/init_visu.jl");
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
cols = cgrad([:black,:blue,:green,:orange,:red,:black]);
plot(rand(10))


L = 200
model1 = VisionXY(L,0.1,pi,"polar")
lattice = TriangularLattice(L)
init_thetas!(model1,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])
# init_thetas!(model1,lattice,init="isolated",q=1,type="source")
    plot_theta(model1,lattice)
