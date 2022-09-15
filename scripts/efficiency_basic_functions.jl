using BenchmarkTools,SpecialFunctions,Distributions,StatsPlots

## Basic Functions
x = 225.3
@btime rand() # 5 ns
@btime randn() # 6 ns
@btime exp($x) # 6 ns
@btime log($x) # 7 ns
@btime cos($x) # 5.5 ns
@btime sin($x) # 5.2 ns
@btime besseli(0,$x) # 5.2 ns

## Some Distributions
@btime rand(VonMises()) # 290 ns
@btime rand(Cosine()) # 1000 ns
@btime rand(Cauchy()) # 33 ns
@btime mod(rand(Cauchy()) -pi , 2pi) # 40 ns
@btime mod(2.0*randn() -pi , 2pi) # 11 ns

data = rand(VonMises(),Int(1E5))
    histogram(data,normalize=true)
    plot!(VonMises())

A = 0.9
    data = mod.(rand(Cauchy(-pi,1/A),Int(1E5)),2pi)
    histogram(data,normalize=true)
    plot!(VonMises(pi,A))

## Trying to numerically sample from VonMises efficiently
dx = 1E-3
xx = -0 : dx : 2π
mu = pi ; k = 10
denom = (2π*besseli(0,k))
vm(x::Float64) = exp(k*cos(x-mu))/denom
yy = f.(xx;mu=π,k=1)
plot(xx,yy)
using QuadGK
using Trapz
start = xx[1]
@btime cdf_vm = [quadgk(vm,$start,x,rtol = dx )[1] for x in xx]
plot(xx,cdf_vm)
u = rand()
@btime xx[findfirst(x->x>u,cdf_vm)]
# cette technique est 100x plus lente
