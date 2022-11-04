include("../src/models.jl");

# Physical Parameters
include(srcdir("../parameters.jl"));

model = XY(params)
model = MonteCarloXY(params)
model = ForcedXY(params)
model = VisionXY(params)
model = SPP(params)

# WARNING, incohérence !
iszero(Float32(π) - (π)) # true
Float32(π) == π # false
