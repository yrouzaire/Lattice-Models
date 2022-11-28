# To scan over
Ts = [0.1,0.2,0.3]
visions = collect(0:0.05:0.2)
scanned_params = Dict("Ts"=>Ts,"visions"=>visions)

# Physical Parameters
L = 50
T = 0.1
symmetry = "polar"
propulsion = "polar"
Var = 0.
A = 0
vision = 0
rho = 1
algo = "A" # rule for collision!() for model = SPP


# Numerical Parameters
dt = 1E-2
float_type = Float32
tmax = 1E2
transients = 1E3
every = 1E2

# Initialisation
init = "hightemp"
q = 1
r0 = Int(L/2)
type1defect = "clockwise"
type2defect = "pair2"

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)
