# LatticeModels

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> LatticeModels

It is authored by Ylann Rouzaire.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

------------------------------------------------------------------------------
# General philosophy of this package
This code aims at simulating the temporal evolution of a 2D vector-field,
parametrized in space by a scalar $\theta(x,y) \in [0,2\pi[$.

Throughout the entire project, we discretize this $\theta$ field onto a discrete 2D lattice.
Pros: the numerical simulations are much faster and easily controlled, in particular
the topological defect detection/tracking.
Cons: discrete lattices tend to induce undesired discretization artifacts, especially the
usual `SquareLattice`.



### Summary:
Each physical model will thus be instantiated and then run following this procedure:
  - Modify parameters.jl with the desired values and load those parameters through the command include(srcdir("../parameters.jl"));.
  - Declare your lattice: SquareLattice(L) or TriangularLattice(L).
  - Declare your model: XY(params), VisionXY(params), SPP(params) etc
  - Initialisation of the $theta$ field: thetas = init_thetas(model,lattice,params_init=params_init)
  - Update the model:
    - For one time step only: update!(thetas,model,lattice)
    - For a given duration $\Delta t$: update!(thetas,model,lattice,$\Delta t$)
    - Until a given time $t$: update!(thetas,model,lattice,tmax=$t$)

# Important files to get started
Here we comment the most important methods in each of these files.
## parameters.jl
Contains physical, numerical and initialisation parameters.
### Physical Parameters
  - $L$ is the (integer) size of the 2D $L\times L$ lattice. It thus contains $L^2$ agents/spins/particles/rotors (all synonyms),
  such that the runtime is of minimum complexity $_mathcal{O}(L^2)$.
  - $T\ge 0$ is the temperature of the bath related to the angle diffusion.
  - `symmetry` is the symmetry of the interaction and can be either "nematic" or "polar".
  Polar means that 2 interacting particles $i$ and $j$ want to point in the same direction, and the force $\propto \sin(\theta_j - \theta_i)$
  Nematic means that 2 interacting particles $i$ and $j$ want to point in the same direction but without the notion of head/tail, and the force $\propto \sin(2(\theta_j - \theta_i))$
  - `propulsion` is the symmetry of propulsion of Self Propelled Particles (SPP model) and can only be "polar" for now. Polar propulsion: self-propelled rods. Nematic propulsion (reverses from times to times): active nematics
  - `Var` $\ge 0$ is the variance of the distribution of intrinsic/internal frequency for the model ForcedXY.
  - $A\ge 0$ is the value of the coupling between orientation and propulsion. $A = 0$: no coupling, $A \to \infty$: high coupling, i.e you are sure you move in the direction of your arrow (if polar propulsion)
  - $0 <\rho \le 1$ is the density. If $\rho 0.9$, 10% of the sites will be empty.
  - `algo` is used independently (but never at the same time) for two cases.
    - If the model is of type XY, `algo = "Langevin"` will make its temporal evolution follow a Langevin dynamics.
    - If the model is of type XY, `algo = "MonteCarlo"` (alias `MC`)  will make its temporal evolution follow a Metropolis MonteCarlo dynamics.
    - If the model is `SPP`, set `algo = "A"` (only one functional)

### Numerical Parameters
  - `dt` is the timestep used in the evolution of `ForcedXY` and `LangevinXY` models.
  - `float_type = Float32` is the recommanded Float type. Float64 takes more memory. Float16 is not a native Float type so no operations are usually slower.
  - `tmax` is the maximum time for temporal evolution.
  - `transients` is the transient time, before which measurements on the system are not conducted.
  - `every` is the interval at which measurements are made, in the case you want linearly spaced data points.

### Initialisation Parameters
  - `init` is the type of initialisation. Can be:
    - "ordered" or "lowtemp" for an ordered initial configuration ($T=0$).
    - "disordered" or "hightemp" for a disordered initial configuration ($T=\infty$).
    - "single" or "isolated" for a manually created vortex (in this case, the charge $q$ and its type `type1defect` should be provided)
    - "pair" for a manually created vortex pair (in this case, the absolute value of the charge $q$, the interdefect separation `r0` (integer) and its type `type2defect` should be provided)
    - "2pair" or "2pairs" for TWO manually created vortex pairs (in cross shape).

## lattice.jl

## models.jl
## init_visu.jl
## core_methods.jl

# Other files
## auxiliary.jl
## measurements.jl

# About topological defects (detection, tracking in-vivo)
The code is located in the defects_methods.jl file.
I will for now only comment on the first 3 methods, used to detect defects.
