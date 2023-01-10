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
r0 = 32#round(Int,L/4)
mu0 = pi/8 # for one defect only
mu_plus,mu_minus,phi = pi/2,nothing,0 # one of them has to be Nothing

# Containers
params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"algo"=>algo)
params_num  = Dict("dt"=>dt,"float_type"=>float_type,"tmax"=>tmax,"transients"=>transients)
params_init = Dict("init"=>init,"q"=>q,"r0"=>r0,"phi"=>phi,"mu0"=>mu0,"mu_minus"=>mu_minus,"mu_plus"=>mu_plus)

params = merge(params_num,params_phys,params_init)
