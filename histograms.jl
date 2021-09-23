include("requirements.jl")
NSIM = 10^4
p0 = parameters()
function HabitatLayer(;height = 1000,color = colorant"red")
    layer(xmin = [0],xmax = [1],ymin = [0],ymax = [height],
    color = [color],
    Geom.rect)
end

swims = [
    (vm_generic,[1.,1//4]),
    (vm_generic,[5//12,1.]),
    (vm_generic,[5//12,1//4]),
    p -> vm_neutral(p,13//24,1//24)
]
swim_names = [:ovm,:dvm,:hybrid,:passive]

V = DataFrame(
    :behavior => repeat(swim_names,inner = NSIM),
    :Xf => zeros(4*NSIM),
    :Xf2 => zeros(4*NSIM),
    :fate => fill(:alex,4*NSIM)
)
H = copy(V)

Threads.@threads for m in 1:4

    if m < 4
        z = swims[m][1](p0,swims[m][2]...)
    else
        z = swims[4]
    end

    L = larva(p0,fv,z,NSIM)
    V.Xf[(1:NSIM) .+ (m-1)*NSIM] = L.stats.Xf
    V.Xf2[(1:NSIM) .+ (m-1)*NSIM] = L.stats.Xf
    V.fate[(1:NSIM) .+ (m-1)*NSIM] = L.stats.fate

    L = larva(p0,fh,z,NSIM)
    H.Xf[(1:NSIM) .+ (m-1)*NSIM] = L.stats.Xf
    H.Xf2[(1:NSIM) .+ (m-1)*NSIM] = L.stats.Xf
    H.fate[(1:NSIM) .+ (m-1)*NSIM] = L.stats.fate
end

for k = 1:nrow(V)
    if V.Xf2[k] <= 1 && V.fate[k] == :settled
        V.Xf2[k] = rand()
    end
    if H.Xf2[k] <= 1 && V.fate[k] == :settled
        H.Xf2[k] = rand()
    end
end



dfhyb = filter(r -> (r.behavior == :hybrid) & (r.fate != :settled),V)
hyb = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfhyb,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfhyb.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfdvm = filter(r -> (r.behavior == :dvm) & (r.fate != :settled),V)
dvm = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfdvm,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfdvm.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfovm = filter(r -> (r.behavior == :ovm) & (r.fate != :settled),V)
ovm = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfovm,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfovm.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfpassive = filter(r -> (r.behavior == :passive) & (r.fate != :settled),V)
pas = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfpassive,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfpassive.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    # Scale.color_discrete_hue(begin dummy(a) = ifelse(a == :dead,colorant"gray",colorant"black"); end),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

draw(PDF(figdir*"/hist-hybrid-0720-hab.pdf",12.5cm,8cm),hyb)
draw(PDF(figdir*"/hist-ovm-0720-hab.pdf",12.5cm,8cm),ovm)
draw(PDF(figdir*"/hist-dvm-0720-hab.pdf",12.5cm,8cm),dvm)
draw(PDF(figdir*"/hist-passive-0720-hab.pdf",12.5cm,8cm),pas)






dfhyb2 = filter(r -> (r.behavior == :hybrid) & (r.fate != :settled),H)
hyb2 = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfhyb2,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfhyb.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfdvm2 = filter(r -> (r.behavior == :dvm) & (r.fate != :settled),H)
dvm2 = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfdvm2,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfdvm.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfovm2 = filter(r -> (r.behavior == :ovm) & (r.fate != :settled),H)
ovm2 = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfovm2,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfovm.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

dfpassive2 = filter(r -> (r.behavior == :passive) & (r.fate != :settled),H)
pas2 = plot(
    AxisLayers([0,60],[0,500])...,
    HabitatLayer(color = colorant"yellowgreen"),
    layer(dfpassive2,x = :Xf,color = :fate,Geom.histogram(bincount = Int(ceil(maximum(dfpassive.Xf))))),
    Coord.cartesian(xmin = 0,xmax = 60,ymin = 0,ymax = 500),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    # Scale.color_discrete_hue(begin dummy(a) = ifelse(a == :dead,colorant"gray",colorant"black"); end),
    Scale.color_discrete_manual(RGB(.1,.1,.1),RGB(.8,.8,.8),levels = [:dead,:lost]),
    my_theme_big_key
)

draw(PDF(figdir*"/ns-hist-hybrid-0721-hab.pdf",12.5cm,8cm),hyb2)
draw(PDF(figdir*"/ns-hist-ovm-0721-hab.pdf",12.5cm,8cm),ovm2)
draw(PDF(figdir*"/ns-hist-dvm-0721-hab.pdf",12.5cm,8cm),dvm2)
draw(PDF(figdir*"/ns-hist-passive-0721-hab.pdf",12.5cm,8cm),pas2)
