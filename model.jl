mutable struct NDProcess
    f::Function
    g::Function
    μ::Function
end

struct LarvaOutput
    p::Dict
    F::NDProcess
    sim::DataFrame
    stats::DataFrame
    runTime::Float64
end

struct LarvaeOutput
    p::Dict
    F::NDProcess
    sims::DataFrame
    stats::DataFrame
    meta::DataFrame
    runTime::Float64
end

function larva(
    p::Dict,
    F::NDProcess,
    z::Array{Float64,1};
    x0::Float64 = rand(),
    dW::Array{Float64,1} = randn(p[:NT]-1)*sqrt(p[:dt]),
    )

    t0 = time()

    sim = DataFrame(
        :t => p[:tvec],
        :x => x0*ones(p[:NT]),
        :z => z
    )

    flag = 0 # 0 = loss
    final_n = p[:NT] # unless something else happens...
    for n = 1:(p[:NT]-1)

        t = sim.t[n]
        X = sim.x[n]
        Z = sim.z[n]

        # do we settle? do we die?
        if t >= p[:Tc] && X <= 1
            flag = 1 # 1 = settle
            final_n = n+1
            break
        else
            u = rand()
            p_death = 1 - exp(-p[:dt]*F.μ(p,t,X,Z))
            if u < p_death
                flag = 2 # 2 = death
                final_n = n+1
                break
            end
        end

        # otherwise, simulate movement in the next step
        sim.x[n+1] = abs(X + F.f(p,t,X,Z)*p[:dt] + F.g(p,t,X,Z)*dW[n])

    end

    tf = float(sim.t[final_n])
    time_surface = float(sum(z[1:final_n])*p[:dt])
    time_offshore = float(count(sim.x[1:final_n] .>= 1.)*p[:dt])
    time_offshore5 = float(count(sim.x[1:final_n] .>= 5.)*p[:dt])
    N_up = count(diff(z[1:final_n]) .== 1.)
    N_down = count(diff(z[1:final_n]) .== -1.)
    sim[(final_n):p[:NT],:] .= NaN

    stats = DataFrame(
        :fate => ifelse(flag > 0,ifelse(flag == 1,:settled,:dead),:lost),
        :tf => tf,
        :Xf => sim.x[final_n-1],
        :time_surface => time_surface,
        :time_bottom => tf - time_surface,
        :time_offshore => time_offshore,
        :time_offshore5 => time_offshore5,
        :time_nearshore => tf - time_offshore,
        :N_up => N_up,
        :N_down => N_down,
        :N => N_up + N_down,
        :final_n => final_n
    )

    return(LarvaOutput(p,F,sim,stats,time() - t0))

end

function larva(
    p::Dict,
    F::NDProcess,
    zfun::Function;
    dW::Array{Float64,1} = randn(p[:NT]-1)*sqrt(p[:dt]),
    )

    return(larva(p,F,zfun(p),x0 = x0,dW = dW))

end

function larva(
    p::Dict,
    F::NDProcess,
    z::Array{Float64,1},
    NSIM::Int64;
    DW::Array{Float64,2} = randn(p[:NT]-1,NSIM)*sqrt(p[:dt])
    )

    t0 = time()

    # initialize with single simulation
    L0 = larva(p,F,z,dW = DW[:,1])

    sims = repeat(L0.sim,outer = NSIM)
    sims[!,:n] = repeat(1:NSIM,inner = p[:NT])
    stats = repeat(L0.stats,NSIM)
    stats[!,:n] = 1:NSIM


    # sims = DataFrame([name => Vector{typeof(L0.sim[1,name])}(undef,NSIM*p[:NT]) for name in Symbol.(names(L0.sim))])
    # # sims = DataFrame(eltype.(eachcol(L0.sim)),names(L0.sim),NSIM*p[:NT])
    # sims[1:p[:NT],names(L0.sim)] = L0.sim
    # sims[!,:n] = repeat(1:NSIM,inner = p[:NT])

    # stats = DataFrame([name => Vector{typeof(L0.stats[1,name])}(undef,NSIM) for name in Symbol.(names(L0.stats))])
    # # stats = DataFrame(eltype.(eachcol(L0.stats)),names(L0.stats),NSIM)
    # stats[1,Symbol.(names(L0.stats))] = L0.stats
    # stats[!,:n] = 1:NSIM

    # run lots of simulations
    # Threads.@threads for n = 2:NSIM
    for n = 2:NSIM
        L = larva(p,F,z,dW = DW[:,n])
        sims[(1:p[:NT]) .+ p[:NT]*(n-1),propertynames(L0.sim)] = L.sim
        stats[n,propertynames(L0.stats)] = L.stats[1,:]
    end

    Ssims = filter(r -> r.fate == :settled,stats)
    meta = DataFrame(
        :settled => count(stats.fate .== :settled)/NSIM,
        :dead => count(stats.fate .== :dead)/NSIM,
        :lost => count(stats.fate .== :lost)/NSIM,
    )
    for metric in names(L0.stats)[2:end]
        name = Symbol("ES_",metric)
        meta[!,name] = [mean(Ssims[:,metric])]
    end

    return(LarvaeOutput(p,F,sims,stats,meta,time() - t0))
end


function larva(
    p::Dict,
    F::NDProcess,
    zfun::Function,
    NSIM::Int64;
    DW::Array{Float64,2} = randn(p[:NT]-1,NSIM)*sqrt(p[:dt])
    )

    t0 = time()

    # initialize with single simulation
    L0 = larva(p,F,zfun(p),dW = DW[:,1])

    sims = repeat(L0.sim,outer = NSIM)
    sims[!,:n] = repeat(1:NSIM,inner = p[:NT])
    stats = repeat(L0.stats,NSIM)
    stats[!,:n] = 1:NSIM


    # sims = DataFrame([name => Vector{typeof(L0.sim[1,name])}(undef,NSIM*p[:NT]) for name in Symbol.(names(L0.sim))])
    # # sims = DataFrame(eltype.(eachcol(L0.sim)),names(L0.sim),NSIM*p[:NT])
    # sims[1:p[:NT],names(L0.sim)] = L0.sim
    # sims[!,:n] = repeat(1:NSIM,inner = p[:NT])

    # stats = DataFrame([name => Vector{typeof(L0.stats[1,name])}(undef,NSIM) for name in Symbol.(names(L0.stats))])
    # # stats = DataFrame(eltype.(eachcol(L0.stats)),names(L0.stats),NSIM)
    # stats[1,Symbol.(names(L0.stats))] = L0.stats
    # stats[!,:n] = 1:NSIM

    # run lots of simulations
    # Threads.@threads for n = 2:NSIM
    for n = 2:NSIM
        L = larva(p,F,zfun(p),dW = DW[:,n])
        sims[(1:p[:NT]) .+ p[:NT]*(n-1),propertynames(L0.sim)] = L.sim
        stats[n,propertynames(L0.stats)] = L.stats[1,:]
    end

    Ssims = filter(r -> r.fate == :settled,stats)
    meta = DataFrame(
        :settled => count(stats.fate .== :settled)/NSIM,
        :dead => count(stats.fate .== :dead)/NSIM,
        :lost => count(stats.fate .== :lost)/NSIM,
    )
    for metric in names(L0.stats)[2:end]
        name = Symbol("ES_",metric)
        meta[!,name] = [mean(Ssims[:,metric])]
    end

    return(LarvaeOutput(p,F,sims,stats,meta,time() - t0))
end
