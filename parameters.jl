# Physical Parameters
L = 200
T = 0.1
symmetry = "polar"
propulsion = "polar"
Var = 0.1
A = 1.
vision = 4π/3
rho = 0.99
antiferro = true


# Numerical Parameters
dt = 1E-2
float_type = Float32
width_proposal = 0.01

# Initialisation
init = "single"
q = 1
r0 = Int(L/2)
type1defect = "clockwise"
type2defect = ["source" , "convergent"]

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"antiferro"=>antiferro)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"width_proposal"=>width_proposal)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)
