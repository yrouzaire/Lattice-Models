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
    params["q"] = 1/2
    params["rho"] = 1
    params["type1defect"] = "counterclockwise"
    params["symmetry"] = "polar"
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    p = plot_thetas(thetas,model,lattice,defects=true)
    window = 9
    display_quiver!(p,thetas,window)

## Plot the different defects and defect pairs
include(srcdir("../parameters.jl"));
    cols = cgrad([:black,:blue,:green,:orange,:red,:black])
    params["L"] = 20
    windows = Int(params["L"]/2-1)
    model = XY(params) # in fact, useless for plotting defects at t = 0
    lattice = SquareLattice(L,periodic=true,single=true) # in fact, useless for plotting defects at t = 0

# +1 Defects
plotsP1 = []
    params["q"] = +1
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["source","sink","clockwise","counterclockwise"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="+1 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsP1,p)
    end
    pP1 = plot(plotsP1...,layout=(1,4),size=(400*4,400))

# -1 Defects
plotsM1 = []
    params["q"] = -1
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["convergent","divergent","threefold1","threefold2"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="-1",size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsM1,p)
    end
    pM1 = plot(plotsM1...,layout=(1,4),size=(400*4,400))

# +1/2 Defects
plotsP12 = []
    params["q"] = +1/2
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["source","sink","clockwise","counterclockwise"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="+1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsP12,p)
    end
    pP12 = plot(plotsP12...,layout=(1,4),size=(400*4,400))

# -1/2 Defects
plotsM12 = []
    params["q"] = -1/2
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["join","split","threefold1","threefold2"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="-1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsM12,p)
    end
    pM12 = plot(plotsM12...,layout=(1,4),size=(400*4,400))

# 4 ≠ Pairs (all the others are equivalent to one of those)
L = 20 ; init = "pair" ; rho = 1
    r0 = 10 ; T = 1 # dummies
    plotsPairs = []
    defect_plus_type = "sink"
    for type in [[defect_plus_type,"convergent"],[defect_plus_type,"divergent"],[defect_plus_type,"threefold1"],[defect_plus_type,"threefold2"]]
        q = 1/2
        thetas = initialize(L,init,q,type,r0,rho)
        p1 = heatmap(mod.(thetas,π),c=cols,clims=(0,π),size=(512,512),aspect_ratio=1,axis=false)
        display_quiver(p1,thetas,window) ; xlims!(1,2window+1) ; ylims!(1,2window+1)
        title!(type[1]*" & "*type[2])

        p2 = heatmap(mod.(thetas,2π),c=cols,clims=(0,2π),size=(512,512),colorbar=false,aspect_ratio=1,axis=false)
        display_quiver(p2,thetas,window) ; xlims!(1,2window+1) ; ylims!(1,2window+1)

        p=plot(p1,p2,layout=(2,1),size=(1024,512))
        push!(plotsPairs,p)
    end
pPairs2 = plot(plotsPairs...,layout=(1,4),size=(1600,800))

# savefig(pP1,"figures/illustration_defects/P1_defects.png")
# savefig(pP1,"figures/illustration_defects/P1_defects.svg")
#
# savefig(pP12,"figures/illustration_defects/P12_defects.png")
# savefig(pP12,"figures/illustration_defects/P12_defects.svg")
#
# savefig(pM1,"figures/illustration_defects/M1_defects.png")
# savefig(pM1,"figures/illustration_defects/M1_defects.svg")
#
# savefig(pM12,"figures/illustration_defects/M12_defects.png")
# savefig(pM12,"figures/illustration_defects/M12_defects.svg")
#
# savefig(pPairs1,"figures/illustration_defects/Pairs1_defects.png")
# savefig(pPairs1,"figures/illustration_defects/Pairs1_defects.svg")
#
# savefig(pPairs2,"figures/illustration_defects/Pairs2_defects.png")
# savefig(pPairs2,"figures/illustration_defects/Pairs2_defects.svg")
#
# savefig(pPairs3,"figures/illustration_defects/Pairs3_defects.png")
# savefig(pPairs3,"figures/illustration_defects/Pairs3_defects.svg")
#
# savefig(pPairs4,"figures/illustration_defects/Pairs4_defects.png")
# savefig(pPairs4,"figures/illustration_defects/Pairs4_defects.svg")
