# Physical Parameters
L = 200
T = 0.1
symmetry = "polar"
propulsion = "polar"
Var = 0.
A = 0
vision = 0.1
rho = 1
algo = "Langevin" # rule for collision!() for model = SPP : algo = "A", or type of XY model : algo = "MonteCarlo"/"MC" or"Langevin

# Numerical Parameters
dt = 1E-2
float_type = Float32
tmax = 1E3
transients = 1E3
every = 1E2

# Initialisation
init = "hightemp"
q = -1/2
r0 = round(Int,L/2)
type1defect = "counterclockwise"
type2defect = "pair1"

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)
