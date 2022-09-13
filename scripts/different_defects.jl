using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

## Plot the different defects and defect pairs
include(srcdir("../parameters.jl"));
    cols = cgrad([:black,:blue,:green,:orange,:red,:black])
    params["L"] = 20
    window = Int(params["L"]/2-1)
    model = XY(params) # in fact, useless for plotting defects at t = 0
    lattice = SquareLattice(L,periodic=true,single=true) # in fact, useless for plotting defects at t = 0

# +1 Defects
plotsP1 = []
    params["q"] = +1
    params["symmetry"] = "polar"
    params["init"] = "single"
    types = ["source","sink","clockwise","counterclockwise"]
    for type in types
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="+1 "*type)
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
plots_pairs = []
    params["r0"] = 8
    params["q"] = 1/2
    params["symmetry"] = "polar"
    params["init"] = "pair"
    for type in ["pair1","pair2","pair3","pair4"]
        params["type2defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plots_pairs,p)
    end
    pPairs = plot(plots_pairs...,layout=(1,4),size=(1600,400))

# savefig(pP1,plotsdir("illustration_defects/P1_defects.png"))
# savefig(pP1,plotsdir("illustration_defects/P1_defects.svg"))
#
# savefig(pP12,plotsdir("illustration_defects/P12_defects.png"))
# savefig(pP12,plotsdir("illustration_defects/P12_defects.svg"))
#
# savefig(pM1,plotsdir("illustration_defects/M1_defects.png"))
# savefig(pM1,plotsdir("illustration_defects/M1_defects.svg"))
#
# savefig(pM12,plotsdir("illustration_defects/M12_defects.png"))
# savefig(pM12,plotsdir("illustration_defects/M12_defects.svg"))
#
# savefig(pPairs,plotsdir("illustration_defects/Pairs_defects.png"))
# savefig(pPairs,plotsdir("illustration_defects/Pairs_defects.svg"))

## Discriminate between same q but different defects : divergence and rotationnal ?
# First implement, then test on my defects, then export to recognition by DefectTracker
include(srcdir("../parameters.jl"));
lattice = SquareLattice(L)
thetas = init_thetas(lattice,params=params)
model = XY(params)
window = 9 # for L = 20


p = plot_thetas(thetas,model,lattice)
    display_quiver!(p,thetas,window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)

divergence,rotational = get_div_rot(thetas)
heatmap(divergence',aspect_ratio=1,size=(485,400))
heatmap(rotational',aspect_ratio=1,size=(485,400))
p

get_divergence(get_neighbours(thetas,model,lattice,10,10))
get_rotational(get_neighbours(thetas,model,lattice,10,10))
