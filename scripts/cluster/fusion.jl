using JLD,Parameters

R = 3
base_filename = "data/PhaseSpace_L$(L)_tmax$(tmax)"
data = load(base_filename*"_r1.jld")
@unpack phi,As,Ts,L,dt,tmax,tsave = data

S = zeros(length(As),length(Ts),nsave,R)
C = zeros(length(As),length(Ts),Lover2,nsave,R)

for r in 1:R
    S[:,:,:,r],Ctw[:,:,:,:,r] = load(base_filename*"_r$r.jld","S","C")
end

# copy all the data (T,Var,L etc) from a random file
filename_fusion = base_filename*"_r1.jld"
# load this base file
data = load(filename_fusion)
# add it all the extra info
data["C"] = C
data["S"] = S
data["R"] = R
# save the whole thing
JLD.save(filename_fusion,data)
println("Fusionned data saved in $filename_fusion.")
