include("../src/models.jl");

# Physical Parameters
include(srcdir("../parameters.jl"));

model = XY(params_phys,params_num)
model = AXY(params_phys,params_num)
model = VisionXY(params_phys,params_num)
model = MovingXY(params_phys,params_num)

# WARNING, incohérence !
iszero(Float32(π) - (π)) # true
Float32(π) == π # false
