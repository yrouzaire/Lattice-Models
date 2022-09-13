# Physical Parameters
L = 200
T = 0.1
symmetry = "nematic"
propulsion = "polar"
Var = 0.2
A = 1.
vision = 5π/3
rho = 1
algo = "A" # rule for collision!() for model = MovingXY

# Numerical Parameters
dt = 1E-2
float_type = Float32
width_proposal = 0.01

# Initialisation
init = "pair"
q = 1
r0 = Int(L/2)
type1defect = "source"
type2defect = "pair1"

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"width_proposal"=>width_proposal)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)
