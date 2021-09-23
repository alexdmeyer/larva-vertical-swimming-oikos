include("requirements.jl")

N = 16
NSIM = 10000
swims = [
    (vm_generic,[1.,1//4]),
    (vm_generic,[5//12,1.]),
    (vm_generic,[5//12,1//4]),
    p -> vm_neutral(p,13//24,1//24)
]
swim_names = [:ovm,:dvm,:hybrid,:passive]

dfu = DataFrame(
    :scheme => repeat([:v,:h],inner = 4N),
    :behavior => repeat(swim_names,inner = N,outer = 2),
    :u1 => repeat(range(0,10,length = N),outer = 8),
    :settled => zeros(8N),
    :settled_h => zeros(8N),
)

dfk = DataFrame(
    :scheme => repeat([:v,:h],inner = 4N),
    :behavior => repeat(swim_names,inner = N,outer = 2),
    :k1 => repeat(10 .^ range(0,2,length = N),outer = 8),
    :settled => zeros(8N),
)

dfT = DataFrame(
    :scheme => repeat([:v,:h],inner = 4N),
    :behavior => repeat(swim_names,inner = N,outer = 2),
    :T => repeat(10 .^ range(0,2,length = N),outer = 8),
    :settled => zeros(8N),
)

p0 = parameters()
Threads.@threads for m = 1:4
    if m < 4
        z = swims[m][1](p0,swims[m][2]...)
    else
        z = swims[4]
    end
    for n = 1:N
        # p = parameters(u1 = dfu.u1[n])
        # dfu.settled[n + N*(m-1)] = larva(p,fv,z,NSIM).meta.settled[1]
        # dfu.settled[n + N*(m+3)] = larva(p,fh,z,NSIM).meta.settled[1]
        # p = parameters(k1 = dfk.k1[n])
        # dfk.settled[n + N*(m-1)] = larva(p,fv,z,NSIM).meta.settled[1]
        # dfk.settled[n + N*(m+3)] = larva(p,fh,z,NSIM).meta.settled[1]
        p = parameters(T = dfT.T[n])
        if m < 4
            zT = swims[m][1](p,swims[m][2]...)
        else
            zT = swims[4]
        end
        dfT.settled[n + N*(m-1)] = larva(p,fv,zT,NSIM).meta.settled[1]
        dfT.settled[n + N*(m+3)] = larva(p,fh,zT,NSIM).meta.settled[1]
    end
end

f1 = plot(
    filter(r -> r.scheme == :v,dfu),
    AxisLayers([0,10],[0,1])...,
    x = :u1,
    y = :settled,
    color = :behavior,
    Geom.line,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    # Geom.SubplotGrid(Geom.line),
    my_theme_big_key,
)

draw(PDF(figdir*"/sens-u-curve-0720.pdf",12.5cm,8cm),f1)

f2 = plot(
    filter(r -> r.scheme == :v,dfk),
    AxisLayers([1,100],[0,1])...,
    x = :k1,
    y = :settled,
    color = :behavior,
    Geom.line,
    Scale.x_log10,
    Guide.xticks(ticks = [0,1,2]),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme_big_key,    # Geom.SubplotGrid(Geom.line)
)

draw(PDF(figdir*"/sens-k-curve-0720.pdf",12.5cm,8cm),f2)

f3 = plot(
    filter(r -> r.scheme == :v,dfT),
    AxisLayers([1,100],[0,1])...,
    x = :T,
    y = :settled,
    color = :behavior,
    Geom.line,
    Scale.x_log10,
    Guide.xticks(ticks = [0,1,2]),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme_big_key,    # Geom.SubplotGrid(Geom.line)
)

draw(PDF(figdir*"/sens-k-curve-0720.pdf",12.5cm,8cm),f2)


f1h = plot(
    filter(r -> r.scheme == :h,dfu),
    AxisLayers([0,10],[0,1])...,
    x = :u1,
    y = :settled,
    color = :behavior,
    Geom.line,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    # Geom.SubplotGrid(Geom.line),
    my_theme_big_key,
)

draw(PDF(figdir*"/ns-sens-u-curve-0720.pdf",12.5cm,8cm),f1h)

f2h = plot(
    filter(r -> r.scheme == :h,dfk),
    AxisLayers([1,100],[0,1])...,
    x = :k1,
    y = :settled,
    color = :behavior,
    Geom.line,
    Scale.x_log10,
    Guide.xticks(ticks = [0,1,2]),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme_big_key,    # Geom.SubplotGrid(Geom.line)
)

draw(PDF(figdir*"/ns-sens-k-curve-0720.pdf",12.5cm,8cm),f2h)
