cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

rand(RGB)
RGB(1,1,1)
convert_RGB(x::Vector{<:Number}) = RGB(x[1],x[2],x[3])
invconvert_RGB(x::Vector{<:Number}) = RGB(1-x[1],1-x[2],1-x[3])
function aliment_sources!(ink,locations,radius;colors)
    for n in each(locations)
        x0,y0 = locations[n]
        aliment_source!(ink,x0,y0,radius,color=colors[n])
    end
end
function aliment_source!(ink,x0,y0,radius;color)
    L = size(ink,1)
    for j in 1:L , i in 1:L
        if (i-x0)^2 + (j-y0)^2  ≤ radius^2
            ink[i,j] = color
        end
    end
end
function visualise_information_flow(thetas,lattice,vision,initial_drops,color_drops,radius,fmax,tmax,every,source_lifetime)
    L = size(thetas,1)
    influenced_by = Matrix{Tuple{Int,Int}}(undef,L,L)
    for j in 1:L , i in 1:L
        a,b = project_angle_onto_lattice(thetas[i,j],i,j,lattice)
        influenced_by[i,j] = Int(mod1(i+a,L)) , Int(mod1(j+b,L))
    end
    ink = fill([0.,0.,0.],L,L)
    aliment_sources!(ink,initial_drops,radius,colors=color_drops) # First ink drops

    # Diffusion in time
    nnn = number_nearest_neighbours(lattice)
    flow = fmax/nnn
    anim = @animate for tt in 0:every:tmax
        for e in 1:every
            ink_new = copy(ink)
            for j in 1:L, i in 1:L
                C = ink[i,j] # concentration
                ink_new[i,j] -= nnn*flow*C
                offs = offsets(lattice,iseven(i))
                for off in offs
                    ii,jj = mod1(i+off[1],L) , mod1(j+off[2],L)
                    ink_new[ii,jj] += (1-vision)*flow*C
                end
                ink_new[influenced_by[i,j]...] += (nnn*vision)*flow*C
            end
            if tt ≤ source_lifetime
                aliment_sources!(ink_new,initial_drops,radius,colors=color_drops) # First ink drops
            end
            ink = ink_new
        end
        heatmap(invconvert_RGB.(ink)',aspect_ratio=1,size=(512,512))
    end
    mp4(anim,"films/information_flow/information_flow_vision$(vision)_tmax$(tmax)_$(length(initial_drops))source.mp4")
    return anim
end


include(srcdir("../parameters.jl"));
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
nnn = number_nearest_neighbours(lattice)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
update!(thetas,model,lattice,tmax=200)
    plot_thetas(thetas,model,lattice,defects=false)

location_drops = reverse.([(125,150),(50,75),(150,75)])
location_drops = [(125,150),(50,75),(150,75)]
color_drops = [[0,1,1],[1,0,1],[1,1,0]]
location_drops = [(100,100)]
color_drops = [[0,1,1]/10]
tmax = 2000 ; source_lifetime = 1 ; every = 10
    fmax = 1
    radius_drop = 200
    vis = 3vision
    z = @elapsed visualise_information_flow(thetas,lattice,vis,location_drops,color_drops,radius_drop,fmax,tmax,every,source_lifetime)
    prinz(z)

## Another way to visualise the system
npaths = 2000
    locations_agents = [(rand(1:L),rand(1:L)) for n in 1:npaths]
    # couleurs = ColorSchemes.tab20.colors[1:npaths]
    # initialisation
    ink = ones(RGB,L,L)
    for n in 1:npaths
        ink[locations_agents[n]...] = RGB(0,0,0) #couleurs[n]
    end
    heatmap(ink')

anim = @animate for tt in 1:200
        for n in 1:npaths
            curr_pos = locations_agents[n]
            next_pos = curr_pos .+ project_angle_onto_lattice(thetas[curr_pos...],curr_pos...,lattice)
            next_pos = mod1.(next_pos,L)
            locations_agents[n] = next_pos# current pos
            ink[next_pos...] = RGB(0,0,0) #couleurs[n]
        end
        heatmap(ink',aspect_ratio=1,size=(512,512))
    end
heatmap(ink',aspect_ratio=1,size=(512,512))
plot_thetas(thetas,model,lattice,defects=false)

mp4(anim,"films/information_flow/$(npaths)paths_vision$(vision)_tmax$(tmax).mp4",fps=25)
