cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));

using Augmentor,LsqFit
#=
The idea of this document is to fit a nematic theta field
around a defect by three numbers, q, µ, and the angle of the rotation \alpha.
And then drop \alpha.
=#

@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP12.jld2")
model = XY(params)
lattice = SquareLattice(W21)

# Generate rotated configurations
degree = 90#rand(0:10:350)
    ppl = Rotate(degree) |> Resize(W21,W21)
    thetas = augment(base_dataset[:,:,1],ppl) .+ Float32(deg2rad(degree))
    p=plot_thetas(thetas,model,lattice,defects=false)
    display_quiver!(p,thetas,WINDOW)

model_thetas(xx,p) = p[1] + 1/2*atan(((xx[1]-WINDOW)*sin(p[2]) + (xx[2]-WINDOW)*cos(p[2]))/((xx[1]-WINDOW)*cos(p[2])-(xx[2]-WINDOW)*sin(p[2])))
p0 = [0,0]
xy = [[x,y] for x in 1:W21, y in 1:W21]
fit = curve_fit(model_thetas, xy, thetas, p0)

using DataFitting
f(x, y, p1, p2) = @.  p1 * x  +  p2 * y
dom = CartesianDomain(30, 40)
d = fill(0., size(dom));
for i in 1:length(dom[1])
    for j in 1:length(dom[2])
        d[i,j] = f(dom[1][i], dom[2][j], 1.2, 2.4)
    end
end
data = Measures(d + randn(rng, size(d)), 1.)
