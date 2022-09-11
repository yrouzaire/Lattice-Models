#= CheckList to pass
1. Verify to which displacement and which color corrresponds θ = 0,π/2,π,-π/2 on a Square Lattice.
2. Check that the colors of the visual aspect of a +/- defect rotates the good way.
3. Check that the defects are correctly localized
4. Check that the arrows of quiver! indicate the correct direction.
4bis. Check that the holes are correctly localized.
=#


## Step 1 OK
L = 10
    model = MovingXY(0.0,100.0,"polar","polar",0.0,1.0,"A",0.01) # T,A,symmetry,propulsion,t,rho,algo,width_proposal
    lattice = TriangularLattice(L)
i,j = 6,5 ; theta = 4pi/3
    thetas = NaN*zeros(L,L)
    thetas[i,j] = theta
    display(plot_thetas(thetas,model,lattice))
    NN = angle2neighbour(thetas[i,j],i,j,model,lattice)
    thetas[i,j] = NaN
    thetas[NN...] = theta
    display(plot_thetas(thetas,model,lattice))
#= Results for SquareLattice
θ = 0 is going up
θ = pi/2 is going left
θ = pi is going down
θ = -pi/2 is going right

Results for TriangularLattice
θ = 0 is going up
θ = pi/2 is going left
θ = pi is going down
θ = -pi/2 is going right (as for Square)
but now, for 'directions obliques', it depends on the parity of "i"

=#

# Faire la roue des couleurs

## Step 2 OK
include(srcdir("../parameters.jl"));
    params["init"] = "single"
    params["q"] = 1/2
    params["type1defect"] = "sink"
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    display(plot_thetas(thetas,model,lattice,defects=true))


## Step 3 OK
include(srcdir("../parameters.jl"));
    params["init"] = "single"
    params["q"] = 1/2
    params["symmetry"] = "nematic"
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
x0,y0 = round(Int,L/2),round(Int,L/2)
    thetas = Float32.(create_single_defect(L,x0,y0,q=q,type=type1defect)) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
    display(plot_thetas(thetas,model,lattice,defects=true))
x0,y0 = round(Int,L/2),round(Int,L/4)
    thetas = Float32.(create_single_defect(L,x0,y0,q=q,type=type1defect)) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
    display(plot_thetas(thetas,model,lattice,defects=true))
x0,y0 = round(Int,L/4),round(Int,L/2)
    thetas = Float32.(create_single_defect(L,x0,y0,q=q,type=type1defect)) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
    display(plot_thetas(thetas,model,lattice,defects=true))
x0,y0 = round(Int,L/4),round(Int,L/4)
    thetas = Float32.(create_single_defect(L,x0,y0,q=q,type=type1defect)) # in case of something more exotic, recall that the use is : create_single_defect(q,type,L,y0,x0) (x and y swapped)
    display(plot_thetas(thetas,model,lattice,defects=true))

## Step 4 and 4bis OK alhamdoullah
include(srcdir("../parameters.jl"));
    params["init"] = "single"
    params["q"] = -1/2
    params["rho"] = .95
    params["type1defect"] = "convergent"
    params["symmetry"] = "polar"
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    p = plot_thetas(thetas,model,lattice,defects=true)
    window = 15
    display_quiver!(p,thetas,window)
