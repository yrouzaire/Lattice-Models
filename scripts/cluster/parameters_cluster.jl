# Physical Parameters
L = 30
Ts = [0.1,0.2]
symmetry = "polar"
propulsion = "polar"
Var = 0.
As = [0,3]
vision = 5π/3
rhos = [0.975,0.9]
algo = "A" # rule for collision!() for model = SPP

# Numerical Parameters
dt = 1E-2
float_type = Float32
tmax = 2E3
transients = 1E3
every = 1E2

# Initialisation
init = "hightemp"
q = 1
r0 = Int(L/2)
type1defect = "clockwise"
type2defect = "pair2"

# Containers
params_phys = Dict("L"=>L,"Ts"=>Ts,"Var"=>Var,"As"=>As,"rhos"=>rhos,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)
