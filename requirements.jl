using Pkg
using DataFrames, DataStructures, CSV
using Gadfly, Cairo, Fontconfig, ColorSchemes, Colors
using Random, StatsBase, Statistics
using Contour
using Distributed
addprocs()

include("parameters.jl")
include("model.jl")

# replace this with the directory for your own machine
home = "/Users/alexandermeyer/Dropbox/Research/Project02_VM/Code/vm-revisions"
figdir = home * "/figures"
datadir = home * "/data"

# for plotting
my_theme = Theme(
    major_label_font="CMU Serif",
    major_label_font_size=16pt,
    minor_label_font="CMU Serif",
    minor_label_font_size=14pt,
    key_title_font="CMU Serif",
    key_title_font_size=12pt,
    key_label_font="CMU Serif",
    key_label_font_size=10pt,
    key_swatch_color = colorant"black",
    key_swatch_size = 2mm,
    point_shapes = [Shape.square,Shape.diamond,Shape.utriangle],
    grid_line_width = 0mm,
    line_width = .5mm
    )
my_theme_big_key = Theme(
    major_label_font="CMU Serif",
    major_label_font_size=16pt,
    minor_label_font="CMU Serif",
    minor_label_font_size=14pt,
    key_title_font="CMU Serif",
    key_title_font_size=12pt,
    key_label_font="CMU Serif",
    key_label_font_size=15pt,
    key_swatch_color = colorant"black",
    key_swatch_size = 2mm,
    point_shapes = [Shape.square,Shape.diamond,Shape.utriangle],
    grid_line_width = 0mm,
    line_width = .75mm,
    )
my_theme_no_key = Theme(
    major_label_font="CMU Serif",
    major_label_font_size=16pt,
    minor_label_font="CMU Serif",
    minor_label_font_size=14pt,
    point_shapes = [Shape.square,Shape.diamond,Shape.utriangle],
    grid_line_width = 0mm,
    line_width = .75mm,
    key_position = :none
)
my_theme_border = Theme(
    major_label_font="CMU Serif",
    major_label_font_size=16pt,
    minor_label_font="CMU Serif",
    minor_label_font_size=14pt,
    point_shapes = [Shape.square,Shape.diamond,Shape.utriangle],
    grid_line_width = 0mm,
    line_width = .75mm,
    key_position = :none,
    panel_stroke = colorant"black"
)

gadfly_colors = Scale.color_discrete_hue().f(10)
my_color_subset = [1 6 4 5 2 7];
my_colors = gadfly_colors[my_color_subset];
my_colormap = p -> get(ColorSchemes.linear_blue_5_95_c73_n256,p)
scatter_theme = Theme(discrete_highlight_color = c -> nothing,alphas = [.5],point_size = 1pt,line_width = 0pt);``



# axis_layers(X,Y;thickness = .2mm) = [
#     layer(x = X,y = Y[[1,1]],Geom.line,color = [colorant"black"],Theme(line_width = thickness)),
#     layer(x = X[[1,1]],y = Y,Geom.line,color = [colorant"black"],Theme(line_width = thickness))
#     ]

fv = NDProcess(
    (p,t,x,z) -> ifelse(z == 1,p[:u1],p[:u0]),
    (p,t,x,z) -> ifelse(z == 1,p[:k1],p[:k0]),
    (p,t,x,z) -> ifelse((z == 1) && (.25 <= t % 1 < .75),p[:μ_hi],p[:μ_lo])
)

fv_cbl = NDProcess(
    (p,t,x,z) -> ifelse(z == 1,p[:u1],p[:u0]) * min(x/p[:ℓ],1),
    (p,t,x,z) -> ifelse(z == 1,p[:k1],p[:k0]) * sqrt(min(x/p[:ℓ],1)),
    (p,t,x,z) -> ifelse((z == 1) && (.25 <= t % 1 < .75),p[:μ_hi],p[:μ_lo])
)

fh = NDProcess(
    (p,t,x,z) -> ifelse(z == 1,p[:u1],p[:u0]),
    (p,t,x,z) -> ifelse(z == 1,p[:k1],p[:k0]),
    (p,t,x,z) -> ifelse(x < 1.,p[:μ_hi],p[:μ_lo])
)

fh_cbl = NDProcess(
    (p,t,x,z) -> ifelse(z == 1,p[:u1],p[:u0]) * min(x/p[:ℓ],1),
    (p,t,x,z) -> ifelse(z == 1,p[:k1],p[:k0]) * sqrt(min(x/p[:ℓ],1)),
    (p,t,x,z) -> ifelse(x < 1.,p[:μ_hi],p[:μ_lo])
)

function vm_generic(p,a,b)
    V = zeros(p[:NT])
    V[.!(a/2 .<= (p[:tvec] .% 1.) .< (1 - a/2)) .& (p[:tvec] .< b*p[:T])] .= 1
    return(V)
end

function vm_neutral(p,λ1,λ2) # \lambda = mean length of stay in lower layer
    V = zeros(p[:NT])
    ℓ = 0
    t = 0.
    i = 1
    while t < p[:T]
        i += 1
        tnext = t + ifelse(iseven(i),λ1,λ2)*randexp()
        V[t .<= p[:tvec] .< tnext] .= ℓ
        ℓ = 1 - ℓ
        t = tnext
    end
    return(V)
    # return(Int.(V))
end

p0 = parameters(T = 1.)
L0 = larva(p0,fv,vm_generic(p0,0.,0.),1)

function ContourLayers(df::DataFrame,z::Symbol,c::Float64;color = "white",line_width = 1mm)
    # for now this is specific to this application :( no choice of x,y
    # Avec = unique(df.A); LA = length(Avec)
    # Bvec = unique(df.B); LB = length(Bvec)
    Avec = unique(df.a); LA = length(Avec)
    Bvec = unique(df.b); LB = length(Bvec)
    Z = zeros(LA,LB)
    for i in 1:LA
        for j in 1:LB
            ℓ = LB*(i-1)+j
            Z[i,j] = df[ℓ,z]
        end
    end
    CC = lines(contour(Avec,Bvec,Z,c))
    coords = [coordinates(C) for C in CC]
    layers = [layer(x = xy[1],y = xy[2],Geom.line(preserve_order = true),Theme(default_color = color,line_width = line_width)) for xy in coords]
    return(layers)
end

function AxisLayers(X,Y;thickness = .2mm)
    return((
        layer(x = X,y = Y[[1,1]],Geom.line,color = [colorant"black"],Theme(line_width = thickness)),
        layer(x = X[[1,1]],y = Y,Geom.line,color = [colorant"black"],Theme(line_width = thickness))
    ))
end
