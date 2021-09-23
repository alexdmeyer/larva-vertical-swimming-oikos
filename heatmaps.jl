include("requirements.jl")
NSIM = 10000 # 5000
cols = propertynames(L0.meta) # L0 is from "requirements"

examples = DataFrame(:behavior => [:OVM,:DVM,:HYBRID],:a => [1.,5/12,5/12],:b => [1/4,1.,1/4])
ExamplesLayer(;color = colorant"white") = layer(
    examples,
    x = :a,
    y = :b,
    shape = :behavior,
    Geom.point,
    Theme(
        default_color = color,
        point_shapes = [Shape.square,Shape.diamond,Shape.utriangle],
        point_size = 2mm,
        discrete_highlight_color = c -> RGB(0,0,0),
        key_position = :none))


p = parameters()
p_short = parameters(T = 5.)
p_T = parameters(T = 15.)
p_u = parameters(u1 = 4.)
p_k = parameters(k1 = 8.)

LA = 25
LB = 25
avec = range(0.,1.,length = LA)
bvec = range(0.,1.,length = LB)

passive_transport = p -> vm_neutral(p,13//24,1//24)

frontier = DataFrame(:a => avec[2:end],:b => .5*(p[:u1]/p[:u0] .+ sqrt.((p[:u1]/p[:u0])^2 .+ 4 ./ avec[2:end])))
function FrontierLayer(;color = colorant"white")
    layer(
        frontier,
        x = :a,
        y = :b,
        Geom.line,
        Theme(
            line_style = [:dash],
            default_color = color,
            line_width=1mm
        )
    )
end


df = DataFrame(:a => repeat(avec,inner = LB),:b => repeat(bvec,outer = LA));
for col in cols
    df[!,col] .= 0.
end
df0 = copy(df) # short larval duration
df_h = copy(df) # horizontal mortality instead of vertical
df_vcbl = copy(df)
df_hcbl = copy(df)
df_T = copy(df)
df_u = copy(df)
df_k = copy(df)



t0 = time()
k = 0
Threads.@threads for n = 1:LA*LB
# Distributed.@distributed (+) for n = 1:LA*LB

    L = larva(p,fv,vm_generic(p,df.a[n],df.b[n]),NSIM)
    df[n,cols] = L.meta[1,cols]

    L = larva(p,fv_cbl,vm_generic(p,df.a[n],df.b[n]),NSIM)
    df_vcbl[n,cols] = L.meta[1,cols]

    L = larva(p_short,fv,vm_generic(p_short,df.a[n],df.b[n]),NSIM)
    df0[n,cols] = L.meta[1,cols]

    L = larva(p,fh,vm_generic(p,df.a[n],df.b[n]),NSIM)
    df_h[n,cols] = L.meta[1,cols]

    L = larva(p,fh_cbl,vm_generic(p,df.a[n],df.b[n]),NSIM)
    df_hcbl[n,cols] = L.meta[1,cols]

    L = larva(p_T,fv,vm_generic(p_T,df.a[n],df.b[n]),NSIM)
    df_T[n,cols] = L.meta[1,cols]

    L = larva(p_u,fv,vm_generic(p_u,df.a[n],df.b[n]),NSIM)
    df_u[n,cols] = L.meta[1,cols]

    L = larva(p_k,fv,vm_generic(p_k,df.a[n],df.b[n]),NSIM)
    df_k[n,cols] = L.meta[1,cols]

    global k += 1
    println(string(100k/(LA*LB))*"% complete (thread "*string(Threads.threadid())*")")
    # pct_next = Int(floor(100k/(LA*LB)))
    # if pct_next > pct
    #     pct = pct_next
    #     println(string(pct) * "% complete")
    # end
end
t1 = time()
println(t1 - t0)

df[!,:∂T_s] = (df.settled - df_T.settled)/(p[:T] - p_T[:T])
df[!,:∂u_s] = (df.settled - df_u.settled)/(p[:u1] - p_u[:u1])
df[!,:∂k_s] = (df.settled - df_k.settled)/(p[:k1] - p_k[:k1])
df[!,:∂T_d] = (df.dead - df_T.dead)/(p[:T] - p_T[:T])
df[!,:∂u_d] = (df.dead - df_u.dead)/(p[:u1] - p_u[:u1])
df[!,:∂k_d] = (df.dead - df_k.dead)/(p[:k1] - p_k[:k1])
df[!,:∂T_l] = (df.lost - df_T.lost)/(p[:T] - p_T[:T])
df[!,:∂u_l] = (df.lost - df_u.lost)/(p[:u1] - p_u[:u1])
df[!,:∂k_l] = (df.lost - df_k.lost)/(p[:k1] - p_k[:k1])

# df[!,:ΔT_s] = sign.(df.∂T_s)
# df[!,:Δu_s] = sign.(df.∂u_s)
# df[!,:Δk_s] = sign.(df.∂k_s)
# df[!,:ΔT_d] = sign.(df.∂T_d)
# df[!,:Δu_d] = sign.(df.∂u_d)
# df[!,:Δk_d] = sign.(df.∂k_d)
# df[!,:ΔT_l] = sign.(df.∂T_l)
# df[!,:Δu_l] = sign.(df.∂u_l)
# df[!,:Δk_l] = sign.(df.∂k_l)

ptv0 = larva(p_short,fv,passive_transport,NSIM)
ptv = larva(p,fv,passive_transport,NSIM)
pth0 = larva(p_short,fh,passive_transport,NSIM)
pth = larva(p,fh,passive_transport,NSIM)
ptv0_cbl = larva(p_short,fv_cbl,passive_transport,NSIM)
ptv_cbl = larva(p,fv_cbl,passive_transport,NSIM)
pth0_cbl = larva(p_short,fh_cbl,passive_transport,NSIM)
pth_cbl = larva(p,fh_cbl,passive_transport,NSIM)




## SENSITIVITY

sens1key = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    layer(df,x = :a,y = :b,color = :∂u_s,Geom.rectbin),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.diverging_bwr_20_95_c54_n256,p),
        minvalue = -.1,maxvalue = .1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key,
)
draw(PDF(figdir*"/sens-u-heat-0720-key.pdf",12.5cm,9cm),sens1key)

sens2key = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    layer(df,x = :a,y = :b,color = :∂k_s),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.diverging_bwr_20_95_c54_n256,p),
        minvalue = -.1,maxvalue = .1),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key,
)
draw(PDF(figdir*"/sens-k-heat-0720-key.pdf",12.5cm,9cm),sens2key)

## VERTICAL MORTALITY - KEY METRICS

VS0 = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df0,:settled,ptv0.meta[1,:settled],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df0),
        x = :a,y = :b,color = :settled,
        Theme(key_position = :inside)
        ),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap,minvalue = 0,maxvalue = .75),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-vs-short-0720.pdf",12.5cm,9cm),VS0)


VS = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df,:settled,ptv.meta[1,:settled],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :settled),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap,minvalue = 0,maxvalue = .75),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-vs-0720.pdf",12.5cm,9cm),VS)


VT = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df,:ES_tf,ptv.meta[1,:ES_tf],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_tf),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        # colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-v-tstar-0720.pdf",12.5cm,9cm),VT)


VToffshore = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df,:ES_time_offshore5,ptv.meta[1,:ES_time_offshore5],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_time_offshore5),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
    colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        # colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-v-offshore-0720.pdf",12.5cm,9cm),VToffshore)


VTsurface = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df,:ES_time_surface,ptv.meta[1,:ES_time_surface],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_time_surface),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-v-surface-0720.pdf",12.5cm,9cm),VTsurface)


VN = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_N),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.viridis,p),
        minvalue = 0,
        maxvalue = 50),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/perf-v-n-0720.pdf",12.5cm,9cm),VN)


traitspace = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    my_theme_border
)
draw(PDF(figdir*"/traitspace-0720.pdf",10cm,10cm),traitspace)










## HORIZONTAL MORTALITY - KEY METRICS

HS = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df_h,:settled,pth.meta[1,:settled],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df_h),
        x = :a,y = :b,color = :settled),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap,minvalue = 0,maxvalue = .75),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/ns-perf-vs-0720.pdf",12.5cm,9cm),HS)


HT = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df_h,:ES_tf,pth.meta[1,:ES_tf],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df_h),
        x = :a,y = :b,color = :ES_tf),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        # colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/ns-perf-v-tstar-0720.pdf",12.5cm,9cm),HT)


HToffshore = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df_h,:ES_time_offshore5,pth.meta[1,:ES_time_offshore5],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df_h),
        x = :a,y = :b,color = :ES_time_offshore5),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
    colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        # colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/ns-perf-v-offshore-0720.pdf",12.5cm,9cm),HToffshore)


HTsurface = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    ContourLayers(df_h,:ES_time_surface,pth.meta[1,:ES_time_surface],color = colorant"white")...,
    layer(
        filter(r -> r.settled > 0,df_h),
        x = :a,y = :b,color = :ES_time_surface),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_kryw_5_100_c67_n256,p),
        minvalue = 0,
        maxvalue = 25),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/ns-perf-v-surface-0720.pdf",12.5cm,9cm),HTsurface)


HN = plot(
    ExamplesLayer(color = colorant"white"),
    FrontierLayer(color = colorant"black"),
    layer(
        filter(r -> r.settled > 0,df_h),
        x = :a,y = :b,color = :ES_N),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.viridis,p),
        minvalue = 0,
        maxvalue = 50),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(title = ""),
    my_theme_big_key
)
draw(PDF(figdir*"/ns-perf-v-n-0720.pdf",12.5cm,9cm),HN)

## additional shit






plot(
    df_hcbl,
    # filter(r -> r.settled > 0,df),
    x = :a,
    y = :b,
    color = :settled,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap),
    Geom.rectbin
)

plot(
    filter(r -> r.settled > 0,df_hcbl),
    x = :a,
    y = :b,
    color = :ES_tf,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap),
    Geom.rectbin
)

plot(
    filter(r -> r.settled > 0,df_hcbl),
    x = :a,
    y = :b,
    color = :ES_time_offshore,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap),
    Geom.rectbin
)

plot(
    filter(r -> r.settled > 0,df_hcbl),
    x = :a,
    y = :b,
    color = :ES_time_surface,
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap),
    Geom.rectbin
)









VS0 = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df0),
        x = :a,y = :b,color = :settled),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap,minvalue = 0,maxvalue = .75),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-vs-short-0718-key.pdf",12cm,8cm),VS0)


VS = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :settled),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = my_colormap,minvalue = 0,maxvalue = .75),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-vs-0718-key.pdf",12cm,8cm),VS)


VT = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_tf),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 30),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-v-tstar-0718-key.pdf",12cm,8cm),VT)


VToffshore = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_time_offshore),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 30),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-v-offshore-0718-key.pdf",12cm,8cm),VToffshore)


VTsurface = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_time_surface),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(
        colormap = p -> get(ColorSchemes.linear_tritanopic_krjcw_5_98_c46_n256,p),
        minvalue = 0,
        maxvalue = 30),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-v-surface-0718-key.pdf",12cm,8cm),VTsurface)


VN = plot(
    examples_layer_white,
    frontier_layer_light,
    layer(
        filter(r -> r.settled > 0,df),
        x = :a,y = :b,color = :ES_N),
    Coord.cartesian(xmin = 0,xmax = 1,ymin = 0,ymax = 1),
    Scale.color_continuous(colormap = p -> get(ColorSchemes.linear_green_5_95_c69_n256,p)),
    Geom.rectbin,
    Guide.xlabel(nothing),
    Guide.ylabel(nothing),
    # Guide.colorkey(""),
    my_theme
)
draw(PDF(figdir*"/perf-v-n-0718-key.pdf",12cm,8cm),VN)
