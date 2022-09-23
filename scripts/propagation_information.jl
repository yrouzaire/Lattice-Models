cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

function get_AB(lattice::Abstract1DLattice,r=round(Int,lattice.L/2))
    L = lattice.L
    @assert r ≤ L/2 && iseven(r)

    A = round(Int,(L-r)/2)
    B = round(Int,(L+r)/2)
    return [CartesianIndex(A)],[CartesianIndex(B)]
end

mutable struct Perturbation2
    amplitude::Number
    t0::Number
    duration::Number
    remains::Bool
end
Perturbation2(;amplitude=π/2,remains=true,t0=0,duration=1E-3) = Perturbation2(amplitude,t0,duration,remains)
f(P::Perturbation2,t) = P.amplitude*(1-exp(-(t-P.t0)/P.duration))/(1+exp(-(t-P.t0)/P.duration))

function perturbe!(thetas::AbstractArray,pertu::Perturbation2,locationA)
    pertu.remains ? tmax = Inf : tmax = 2pertu.duration
    if 0 ≤ model.t - pertu.t0 ≤ tmax
        for site in locationA
            thetas[site] = f(pertu,model.t)
        end
    end
    return thetas
end

include(srcdir("../parameters.jl"));
lattice = Chain1D(100,periodic=true)
Perturb = Perturbation2(amplitude=pi,remains=true,duration=50)
# plot(t->f(Perturb,t))
spinsA = CartesianIndex.(collect(40:60))

model = PropagationForcedXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)

while model.t < 2500
    perturbe!(thetas,Perturb,spinsA)
    update!(thetas,model,lattice,model.t+50)
    display(plot_thetas(thetas,model,lattice))
    end

&
# function get_AB(lattice::Abstract2DLattice,r,width=1)
#     L = lattice.L
#     @assert r ≤ L/2 && iseven(r)
#     @assert isinteger(width)
#
#     centerA = (round(Int,(L-r)/2),round(Int,L/2))
#     centerB = (round(Int,(L+r)/2),round(Int,L/2))
#     if width > 1
#         if iseven(width) width += 1 end
#         half_window = (width-1)/2
#         A = ????
#         B = ????
#     else
#         A = centerA
#         B = centerB
#     end
#
    # return CartesianIndex(A),CartesianIndex(B)
# end






## Test update 1D lattice
include(srcdir("../parameters.jl"));
lattice = Chain1D(L,periodic=false)
model = PropagationForcedXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice)

update!(thetas,model,lattice,80)
    plot_thetas(thetas,model,lattice)

## Test update 2D lattice
include(srcdir("../parameters.jl"));
lattice = SquareLattice(L)
model = PropagationForcedXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice)

update!(thetas,model,lattice,30)
    plot_thetas(thetas,model,lattice)
