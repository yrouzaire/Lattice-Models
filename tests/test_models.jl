include("../src/models.jl");

# Physical Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    vision = π
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"vision"=>vision,"symmetry"=>symmetry)
# Numerical Parameters
dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

model = XY(params_phys,params_num)
model = AXY(params_phys,params_num)
model = VisionXY(params_phys,params_num)

# WARNING, incohérence !
iszero(Float32(π) - (π)) # true
Float32(π) == π # false
