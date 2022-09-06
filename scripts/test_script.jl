1+1
using LatticeModels
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)

## Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    propulsion = "polar"
    Var = 0.1
    A = 1.
    vision = 4π/3
    rho = 0.99
    antiferro = true
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry,"propulsion"=>propulsion,"antiferro"=>antiferro)

    dt = 1E-2
    float_type = Float32
    width_proposal = 0.01
    params_num  = Dict("dt"=>dt,"float_type"=>float_type,"width_proposal"=>width_proposal)

## Benchmark update
Normal()
model = MovingXY(params_phys,params_num)
lattice = TriangularLattice(L,periodic=true)
thetas = init_thetas(model,lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
update!(thetas,model,lattice)
