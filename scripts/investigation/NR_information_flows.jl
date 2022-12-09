cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));
model = SoftVisionXY(params)
lattice = TriangularLattice(L)

thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
update!(thetas,model,lattice,tmax=10)
plot_thetas(thetas,model,lattice,defects=false)


L = size(thetas,1)
influenced_by = zeros(Int,L,L)
for j in 1:L , i in 1:L
    a,b = project_angle_onto_lattice(thetas[i,j],i,j,lattice)
    number_in[mod1(i+a,L),mod1(j+b,L)] += 1
end

global ink = ones(RGB,L,L)
# First ink drop
x0,y0 = 50,50
fmax = 0.2 # max flow per unit time
ink[x0,y0] = RGB(1,0,0)
heatmap(ink)

ink_new = copy(ink)
for j in 1:L, i in 1:L
    lost = fmax*ink[i,j]
    ink_new[i,j] -= lost

end


ttmax = 5
for tt in 1:ttmax
    for j in 1:L, i in 1:L
        lost_amount = fmax - nnn*fmax*vision*number_in[i,j]
        r,b,g = red(ink[i,j]),green(ink[i,j]),blue(ink[i,j])
        ink[i,j] = RGB(min(1,r+lost_amount/3),min(1,g+lost_amount/3),min(1,b+lost_amount/3))
    end
    ink[x0,y0] = RGB(1,0,0)
    display(heatmap(ink'))
end
