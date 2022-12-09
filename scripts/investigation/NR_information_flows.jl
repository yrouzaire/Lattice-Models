cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

RGB(x) = RGB(x[1],x[2],x[3])
invRGB(x) = RGB(1-x[1],1-x[2],1-x[3])

function aliment_source!(ink,x0,y0,radius;color)
    L = size(ink,1)
    for j in 1:L , i in 1:L
        if (i-x0)^2 + (j-y0)^2  ≤ radius^2
            ink[i,j] = color
        end
    end
end
function visualise_information_flow(thetas,lattice,initial_drops,color_drops,radius,fmax,tmax)
    L = size(thetas,1)
    influenced_by = Matrix{Tuple{Int,Int}}(undef,L,L)
    for j in 1:L , i in 1:L
        a,b = project_angle_onto_lattice(thetas[i,j],i,j,lattice)
        influenced_by[i,j] = Int(mod1(i+a,L)) , Int(mod1(j+b,L))
    end
    ink = fill([0.,0.,0.],L,L)
    # First ink drops
    for n in each(initial_drops)
        x0,y0 = initial_drops[n]
        aliment_source!(ink,x0,y0,radius,color=color_drops[n])
    end

    # Diffusion in time
    nnn = number_nearest_neighbours(lattice)
    flow = fmax/nnn
    anim = @animate for tt in 1:tmax
        ink_new = copy(ink)
        for j in 1:L, i in 1:L
            C = ink[i,j] # concentration
            lost = nnn*flow*C
            ink_new[i,j] -= lost
            offs = offsets(lattice,iseven(i))
            for off in offs
                ii,jj = mod1(i+off[1],L) , mod1(j+off[2],L)
                ink_new[ii,jj] += (1-vision)*flow*C
            end
            ink_new[influenced_by[i,j]...] += (1+(nnn-1)*vision)*flow*C
        end
        ink = ink_new
        heatmap(invRGB.(ink'))
    end
    mp4(anim,"films/information_flow/test.mp4")
    return anim
end


include(srcdir("../parameters.jl"));
model = SoftVisionXY(params)
lattice = SquareLattice(L)
nnn = number_nearest_neighbours(lattice)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
update!(thetas,model,lattice,tmax=100)
plot_thetas(thetas,model,lattice,defects=false)

location_drops = [(50,50)]
color_drops = [[0,1,1]]
tmax = 10
fmax = 0.4
radius_drop = 2
visualise_information_flow(thetas,lattice,position_drops,color_drops,radius_drop,fmax,tmax)


L = size(thetas,1)
influenced_by = Matrix{Tuple{Int,Int}}(undef,L,L)
for j in 1:L , i in 1:L
    a,b = project_angle_onto_lattice(thetas[i,j],i,j,lattice)
    influenced_by[i,j] = Int(mod1(i+a,L)) , Int(mod1(j+b,L))
end

fmax = 0.4/nnn # max flow per unit time
global ink = fill([0.,0.,0.],L,L)
# First ink drop
x0,y0 = 50,50
radius = 2
ink[x0,y0] = [0,1,1]
for j in 1:L , i in 1:L
    if (i-x0)^2 + (j-y0)^2  ≤ radius^2
        ink[i,j] = [0,1,1]
    end
end

heatmap(invRGB.(ink))

ink_new = copy(ink)
for j in 1:L, i in 1:L
    C = ink[i,j] # concentration
    lost = nnn*fmax*C
    ink_new[i,j] -= lost
    offs = offsets(lattice,iseven(i))
    for off in offs
        ii,jj = mod1(i+off[1],L) , mod1(j+off[2],L)
        ink_new[ii,jj] += (1-vision)*fmax*C
    end
    ink_new[influenced_by[i,j]...] += (1+(nnn-1)*vision)*fmax*C
end
    # Aliment the source
    for j in 1:L , i in 1:L
        if (i-x0)^2 + (j-y0)^2  ≤ radius^2
            ink_new[i,j] = [0,1,1]
        end
    end
    ink = ink_new
    heatmap(invRGB.(ink_new))
plot_thetas(thetas,model,lattice,defects=false)



ttmax = 5
for tt in 1:ttmax
    for j in 1:L, i in 1:L
        lost_amount = fmax - nnn*fmax*vision*number_in[i,j]
        r,b,g = red(ink[i,j]),green(ink[i,j]),blue(ink[i,j])
        ink[i,j] = RGB(min(1,r+lost_amount/3),min(1,g+lost_amount/3),min(1,b+lost_amount/3))
    end
    ink[x0,y0] = RGB(1,0,0)
    display(heatmap(ink'))
end
