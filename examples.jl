include("requirements.jl")
NSIM = 10000;
TMAX = 5.
p = parameters(T = TMAX)
function SettlingBlockLayer(;color = colorant"lightgray")
    layer(
        xmin = [p[:Tc]],
        xmax = [p[:T]],
        ymin = [0],
        ymax = [1],
        color = [color],
        Geom.rect
    )
end

swims = [
    (vm_generic,[1.,1//4]),
    (vm_generic,[5//12,1.]),
    (vm_generic,[5//12,1//4]),
    p -> vm_neutral(p,13//24,1//24)
]
swim_names = [:ovm,:dvm,:hybrid,:passive]

sims = Dict()
lost_larvae = Dict()
Threads.@threads for m in 1:4
    if m < 4
        z = swims[m][1](p,swims[m][2]...)
    else
        z = swims[4]
    end
    sims[swim_names[m]] = larva(p,fv,z,NSIM)
    lost_larvae[swim_names[m]] = findall(sims[swim_names[m]].stats.fate .== :lost)
end




ex_passive_z = plot(
    filter(trt -> (trt.n == lost_larvae[:passive][8]) & !isnan(trt.t),sims[:passive].sims),
    AxisLayers([0,TMAX],[0,1])...,
    x = :t,
    y = :z,
    color = [colorant"black"],
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = TMAX,ymin = 0,ymax = 1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    my_theme_no_key
)

ex_ovm_z = plot(
    filter(trt -> (trt.n == lost_larvae[:ovm][1]) & !isnan(trt.t),sims[:ovm].sims),
    AxisLayers([0,TMAX],[0,1])...,
    x = :t,
    y = :z,
    color = [colorant"black"],
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = TMAX,ymin = 0,ymax = 1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    my_theme_no_key
)

ex_dvm_z = plot(
    filter(trt -> (trt.n == lost_larvae[:dvm][1]) & !isnan(trt.t),sims[:dvm].sims),
    AxisLayers([0,TMAX],[0,1])...,
    x = :t,
    y = :z,
    color = [colorant"black"],
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = TMAX,ymin = 0,ymax = 1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    my_theme_no_key
)

ex_hybrid_z = plot(
    filter(trt -> (trt.n == lost_larvae[:hybrid][1]) & !isnan(trt.t),sims[:hybrid].sims),
    AxisLayers([0,TMAX],[0,1])...,
    x = :t,
    y = :z,
    color = [colorant"black"],
    Geom.line,
    Coord.cartesian(xmin = 0,xmax = TMAX,ymin = 0,ymax = 1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    my_theme_no_key
)



ex_passive_x = plot(
    AxisLayers([0,TMAX],[0,40])...,
    SettlingBlockLayer(color = colorant"black"),
    layer(
        sims[:passive].sims,
        x = :t,
        y = :x,
        Geom.histogram2d(xbincount = Int((p[:NT]-1)/2))),
    Scale.color_log10(colormap = p -> get(ColorSchemes.linear_worb_100_25_c53_n256,p)),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    Coord.cartesian(ymax = 30,xmax = TMAX),
    my_theme_no_key
)



ex_ovm_x = plot(
    AxisLayers([0,TMAX],[0,40])...,
    SettlingBlockLayer(color = colorant"black"),
    layer(
        sims[:ovm].sims,
        x = :t,
        y = :x,
        Geom.histogram2d(xbincount = Int((p[:NT]-1)/2))),
    Scale.color_log10(colormap = p -> get(ColorSchemes.linear_worb_100_25_c53_n256,p)),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    Coord.cartesian(ymax = 30,xmax = TMAX),
    my_theme_no_key
)

ex_dvm_x = plot(
    AxisLayers([0,TMAX],[0,40])...,
    SettlingBlockLayer(color = colorant"black"),
    layer(
        sims[:dvm].sims,
        x = :t,
        y = :x,
        Geom.histogram2d(xbincount = Int((p[:NT]-1)/2))),
    Scale.color_log10(colormap = p -> get(ColorSchemes.linear_worb_100_25_c53_n256,p)),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    Coord.cartesian(ymax = 30,xmax = TMAX),
    my_theme_no_key
)

ex_hybrid_x = plot(
    AxisLayers([0,TMAX],[0,40])...,
    SettlingBlockLayer(color = colorant"black"),
    layer(
        sims[:hybrid].sims,
        x = :t,
        y = :x,
        Geom.histogram2d(xbincount = Int((p[:NT]-1)/2))),
    Scale.color_log10(colormap = p -> get(ColorSchemes.linear_worb_100_25_c53_n256,p)),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    Coord.cartesian(ymax = 30,xmax = TMAX),
    my_theme_no_key
)

ex_hybrid_x_key = plot(
    AxisLayers([0,TMAX],[0,40])...,
    SettlingBlockLayer(color = colorant"black"),
    layer(
        sims[:hybrid].sims,
        x = :t,
        y = :x,
        Geom.histogram2d(xbincount = Int((p[:NT]-1)/2))),
    Scale.color_log10(colormap = p -> get(ColorSchemes.linear_worb_100_25_c53_n256,p)),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    Coord.cartesian(ymax = 30,xmax = TMAX),
    my_theme_big_key
)


draw(PDF(figdir * "/ex-passive-z-0720.pdf",8cm,6cm),ex_passive_z)
draw(PDF(figdir * "/ex-passive-x-0720.pdf",8cm,6cm),ex_passive_x)
draw(PDF(figdir * "/ex-ovm-z-0720.pdf",8cm,6cm),ex_ovm_z)
draw(PDF(figdir * "/ex-ovm-x-0720.pdf",8cm,6cm),ex_ovm_x)
draw(PDF(figdir * "/ex-dvm-z-0720.pdf",8cm,6cm),ex_dvm_z)
draw(PDF(figdir * "/ex-dvm-x-0720.pdf",8cm,6cm),ex_dvm_x)
draw(PDF(figdir * "/ex-hybrid-z-0720.pdf",8cm,6cm),ex_hybrid_z)
draw(PDF(figdir * "/ex-hybrid-x-0720.pdf",8cm,6cm),ex_hybrid_x)
draw(PDF(figdir * "/ex-hybrid-x-0720-key.pdf",10cm,6cm),ex_hybrid_x_key)
