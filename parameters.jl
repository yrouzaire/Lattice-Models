# Physical Parameters
L = 100
T = 0.1
symmetry = "polar"
propulsion = "polar"
Var = 0.1
A = 0
vision = 0.2
rho = 1
algo = "Langevin" # rule for collision!() for model = SPP : algo = "A", or type of XY model : algo = "MonteCarlo"/"MC" or"Langevin

# Numerical Parameters
dt = 1E-1
float_type = Float32
tmax = 1E2
transients = 1E3
every = 1E2

# Initialisation
init = "pair"
q = 1
r0 = 40#round(Int,L/4)
phi = 0
type1defect = 3pi/2
type2defect = "pair3"

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"phi"=>phi,"type1defect"=>type1defect,"type2defect"=>type2defect)

params = merge(params_num,params_phys,params_init)

#
# # Initialisation
# init = "pair"
# q = 1
# r0 = round(Int,L/3)
# phi = 0
# mu0 = 3pi/2 # for one defect
# mu1,mu2,phi = "pair3" # for two defects, one of the 3 parameters has to be nothing
#
# # Containers
# params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
# params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
# params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"phi"=>phi,"mu0"=>mu0,"mu1"=>mu1,"mu2"=>mu2,"phi"=>phi)
#
# params = merge(params_num,params_phys,params_init)
