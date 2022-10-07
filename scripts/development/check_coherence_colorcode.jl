cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


#= CheckList to pass
1. Verify to which displacement and which color corrresponds θ = 0,π/2,π,-π/2 on a Square Lattice.
2. Check that the colors of the visual aspect of a +/- defect rotates the good way.
3. Check that the defects are correctly localized
4. Check that the arrows of quiver! indicate the correct direction.
4bis. Check that the holes are correctly localized.
=#

## Step 0 : Comprendre heatmap
test = reshape(1:16,4,4)
test'
heatmap(test)
heatmap(test)
heatmap(test)
heatmap(test')
heatmap(rotr90(test)')

## Step 1 : OK
# noir vers la droite, blue vers le haut, vert vers la gauche, rouge vers le bas
L = 10
    T = 0.1
    symmetry = "polar"
    propulsion = "polar"
    A = 10
    rho = 5/L^2
    algo = "A" # rule for collision!() for model = MovingXY
    float_type = Float32
    init = "lowtemp"
    q = 1/2
    r0 = round(Int,L/2)
    type1defect = "counterclockwise"
    type2defect = "pair1"

    # Containers
    params_phys = Dict("L"=>L,"T"=>T,"A"=>A,"rho"=>rho,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
    params_num  = Dict("float_type"=>float_type)
    params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)
    params = merge(params_num,params_phys,params_init)

model = MovingXY(params) # T,A,symmetry,propulsion,t,rho,algo,width_proposal
lattice = SquareLattice(L)
thetas = 0*init_thetas(model,lattice,params_init=params_init)
thetas .+= pi/2
plot_thetas(thetas,model,lattice)
update!(thetas,model,lattice)
    plot_thetas(thetas,model,lattice)

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

## Step 2 (ca tourne dans le bon sens) OK
include(srcdir("../parameters.jl"));
params["init"] = "single"
    params["symmetry"] = "polar"
    params_init["q"] = -1/2
    params_init["type1defect"] = pi/2
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(model,lattice,params_init=params_init)
    display(plot_thetas(thetas,model,lattice,defects=true))


## Step 3 OK
include(srcdir("../parameters.jl"));
    params["init"] = "single"
    params["q"] = q = 1/2
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
    params_init["init"] = "single"
    params_init["q"] = -1/2
    params["rho"] = 1
    params_init["type1defect"] = pi
    params["symmetry"] = params_init["symmetry"] = "nematic"
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(model,lattice,params_init=params_init)
    p = plot_thetas(thetas,model,lattice,defects=true)
    window = 9
    zoom_quiver(thetas,model,lattice,100,100,window)
