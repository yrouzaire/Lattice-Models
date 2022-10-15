include("../src/models.jl");

# Physical Parameters
include(srcdir("../parameters.jl"));

model = XY(params)
model = MCXY(params)
model = ForcedXY(params)
model = VisionXY(params)
model = MovingXY(params)

# WARNING, incohérence !
iszero(Float32(π) - (π)) # true
Float32(π) == π # false
