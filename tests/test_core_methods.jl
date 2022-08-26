include("../src/core_methods.jl");

## Tests get neighbours
L = 200
model = XY(L,0.1,"polar")
lattice = TriangularLattice(L)
init_thetas!(model,lattice,init="hightemp")
plot_theta(model,lattice)

i = 1; j=1
get_neighbours(model,lattice,i,j)


L = 200
model = VisionXY(L,0.1,pi-0.1,"polar")
lattice = TriangularLattice(L)
init_thetas!(model,lattice,init="hightemp")
plot_theta(model,lattice)

i = 1; j=1
get_neighbours(model,lattice,i,j)
