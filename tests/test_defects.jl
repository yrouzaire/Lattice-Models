include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");
using BenchmarkTools

## Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    vision = π
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"vision"=>vision,"symmetry"=>symmetry)

dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,init="isolated",q=1,type="source")
# thetas = init_thetas(lattice,init="pair",q=1,r0=60,type=["source","divergent"])
spot_defects(thetas,model,lattice)
plot_theta(thetas,model,lattice,defects=true)


##
abstract type AbstractTest{A,B} end
mutable struct X{A,B} <: AbstractTest{A,B}
    p1
    p2
end
function X(A,B)
    return X{A,B}(1,1)
end
x = X(true,true)
typeof(x)
function f1(x::AbstractTest{A,B}) where {A=Bool,B}
    return "A = true ; B = true"
end
function f2(x::X{A,B}) where {A,B}
    return "General"
end
function f2(x::X{true,true}) where {A,B}
    return "true true"
end
function f2(x::X{A,false}) where {A,B}
    return "A false"
end

f2(X(false,false))


f1(x)
