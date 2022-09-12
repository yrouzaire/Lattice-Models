"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("IDrealisation.jl") ;
using JLD,StatsBase,Distributions,LinearAlgebra,Parameters,Statistics,Random
include("methods.jl");

## Scan parameters
L = 200
rhos = [0.9,0.95,0.96,0.97,0.98,0.99,1.0]
Ts = 0.1
As = [0,0.5,1,5]

antiferro  = true
propulsion = "polar"
init = "hightemp" ; q = 1/2 ; r0 = Int(L/2) ; BC = "periodic" ; type = "source"
tmax = Int(1E6) ; times = unique!(round.(Int,2.0 .^range(0,log2(tmax),length=41)))
stdev = 0.15

comments = "Choice of nearest neighbours with Von Mises. The angle proposals are with wrapped gaussians because stdev² = 1% so both methods would lead to the same results. "

n = zeros(Int,length(rhos),length(Ts),length(As),length(times))
C = zeros(round(Int,L/2),length(rhos),length(Ts),length(As),length(times))

z = @elapsed for i in each(rhos) , j in each(Ts) , k in each(As)
    rho = rhos[i] ; T = Ts[j] ; A = As[k]
    thetas = initialize(L,init,q,type,r0,rho)
    t = 0 ; token  = 2
    n[i,j,k,1]   = number_defects(thetas,L,BC)
    C[:,i,j,k,1] = get_C(thetas)
    while t < tmax
        t += 1
        thetas = base_model(thetas,L,T,A,propulsion,BC,stdev=stdev,antiferro=antiferro)
        if t ≥ times[token]
            n[i,j,k,token]   = number_defects(thetas,L,BC)
            C[:,i,j,k,token] = get_C(thetas)
            token = min(token+1,length(times))
        end
    end
end
prinz(z)

base_filename = "data/testXY"
if antiferro filename = base_filename*"_FAF_r$real.jld"
else         filename = base_filename*"_r$real.jld"
end
JLD.save(filename,"L",L,"Ts",Ts,"As",As,"rhos",rhos,"init",init,"propulsion",propulsion,"BC",BC,"tmax",tmax,"times",times,"n",n,"C",C,"q",q,"type",type,"runtime",z,"comments",comments)
