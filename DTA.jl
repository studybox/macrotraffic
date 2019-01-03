module DTA

using POMDPs
using POMDPToolbox
using ParticleFilters

using Distributions
using StaticArrays
using Parameters
using Interpolations
import ParticleFilters: obs_weight
#=
const Vec3 = SVector{3, Float64}
const Vec4 = SVector{4, Float64}
const Vec4Zeros = Vec4([0,0,0,0])
const Vec4Ones = Vec4([1,1,1,1])
=#
COSTURBAN = [0.1325, 0.1325, 0.1325, 0.1312, 0.1300, 0.1288, 0.1275, 0.1260, 0.1241, 0.1217, 0.1187, 0.1152, 0.1115, 0.1081, 0.1051]
importall POMDPs
import Base: rand

export DTAstate, DTAaction, DTAobservation, DTAMDP, DTAPOMDP, mdp, obs_weight, ObsAdaptiveParticleFilter, generate_sr, isterminal, baseline_action

struct DTAstate
    highway::Float64
    urban::Float64
    uassign::Float64
end

struct DTAaction
    assign :: Float64
end

struct DTAobservation
    highway :: Float64
    urban :: Float64
end

@with_kw immutable DTAMDP <: MDP{DTAstate, DTAaction}
    lambda::Float64 = 200
    shortest::Float64 = 0.6
    connected::Float64 = 0.5
    capacities::Float64= 100
    tau::Float64 = 5*60
    omega::Float64 = 40
    discount::Float64 = 0.99
    itpcost::Any = interpolate(COSTURBAN, BSpline(Linear()), OnGrid())
    cost_queue::Float64 = 0.1338
    cost_highway::Float64 = 0.0861
    queuespeed::Float64 = 12
    highwayspeed::Float64 = 30
    maxurbandensity::Float64 = 0.06
    urbanlength::Float64 = 11000
    highwaylength::Float64 = 12600
    maxurbanspeed::Float64 = 15
end
@with_kw immutable DTAPOMDP <: POMDP{DTAstate, DTAaction, DTAobservation}
    mdp ::DTAMDP = DTAMDP()
    obs_std::Float64 = 0.5
end

const DTAProblem = Union{DTAMDP, DTAPOMDP}

mdp(p::DTAPOMDP) = p.mdp
mdp(p::DTAMDP) = p
#=
function generate_sor(pp::DTAProblem, s::DTAstate, a::DTAaction, rng::AbstractRNG)
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    uassigns = s.uassign
    assigns = a.assign
    new_demands = max(rand(rng, Normal(p.lambda, p.lambda*0.1)), 0)

    qhighway = (p.connected*new_demands.*assigns + (1-p.connected)*new_demands.*uassigns)/p.tau
    u = (qhighway - p.capacities)/(p.capacities/p.queuespeed - qhighway/p.highwayspeed)
    new_highways = min.(max.(highways + u*p.tau, 0), p.highwaylength)

    new_passes = max(rand(rng, Normal(p.omega, p.omega*0.05)), 0)
    enter_urbans = p.connected*new_demands.*(1-assigns) + (1-p.connected)*new_demands.*(1-uassigns)+new_passes
    exit_urbans = Green_field(p, urbans)
    new_urbans = max(enter_urbans +  urbans - exit_urbans, 0)
    o = observation(p, s, a, new_highways, new_urbans, u, rng)
    new_uassign = user_equilibrium(p, o)

    sp = DTAstate(new_highways, new_urbans, new_uassign)
    r = reward(p, s, a, sp, u, new_demands)
    return (sp, o, r)
end
=#
function generate_s(pp::DTAPOMDP, s::DTAstate, a::DTAaction, rng::AbstractRNG)
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    uassigns = s.uassign
    assigns = a.assign
    new_demands = max(rand(rng, Normal(p.lambda, p.lambda*0.1)), 0)
    qhighway = (p.connected*new_demands.*assigns + (1-p.connected)*new_demands.*uassigns)/p.tau
    k = qhighway/p.highwayspeed
    #if qhighway/p.highwayspeed >= 0.8*p.capacities/p.queuespeed/p.tau
    #    k = 0.8*p.capacities/p.queuespeed/p.tau
    #end
    u = (qhighway - p.capacities/p.tau)/(p.capacities/p.tau/p.queuespeed - k)
    new_highways = min.(max.(highways + u*p.tau, 0), p.highwaylength)

    new_passes = max(rand(rng, Normal(p.omega, p.omega*0.05)), 0)
    enter_urbans = p.connected*new_demands.*(1-assigns) + (1-p.connected)*new_demands.*(1-uassigns)+new_passes
    exit_urbans = Green_field(p, urbans)
    new_urbans = max(enter_urbans +  urbans - exit_urbans, 0)
    o = observation(p, s, a, new_highways, new_urbans, u, rng)
    new_uassign = user_equilibrium(p, o)

    sp = DTAstate(new_highways, new_urbans, new_uassign)
    #r = reward(p, s, a, sp, u, new_demands)
    return sp
end

function generate_sr(pp::DTAPOMDP, s::DTAstate, a::DTAaction, rng::AbstractRNG)
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    uassigns = s.uassign
    assigns = a.assign
    new_demands = max(rand(rng, Normal(p.lambda, p.lambda*0.1)), 0)

    qhighway = (p.connected*new_demands.*assigns + (1-p.connected)*new_demands.*uassigns)/p.tau
    k = qhighway/p.highwayspeed
    #if qhighway/p.highwayspeed >= 0.8*p.capacities/p.queuespeed/p.tau
    #    k = 0.8*p.capacities/p.queuespeed/p.tau
    #end
    u = (qhighway - p.capacities/p.tau)/(p.capacities/p.tau/p.queuespeed - k)
    #println(u*p.tau)
    new_highways = min.(max.(highways + u*p.tau, 0), p.highwaylength)
    #println(new_highways)
    new_passes = max(rand(rng, Normal(p.omega, p.omega*0.05)), 0)
    enter_urbans = p.connected*new_demands.*(1-assigns) + (1-p.connected)*new_demands.*(1-uassigns)+new_passes
    exit_urbans = Green_field(p, urbans)
    new_urbans = max(enter_urbans +  urbans - exit_urbans, 0)
    o = observation(p, s, a, new_highways, new_urbans, u, rng)
    new_uassign = user_equilibrium(p, o)

    sp = DTAstate(new_highways, new_urbans, new_uassign)
    r = reward(p, s, a, sp, u, new_demands)
    return (sp, r)
end

function Green_field(p, x)
    density = x/p.urbanlength
    v = max(p.maxurbanspeed*(1-density/p.maxurbandensity), 2.5)
    qu = v*density
    return qu*p.tau
end
function user_equilibrium(p::DTAProblem, o::DTAobservation)
    highwayspeed = o.highway
    urbanspeed = o.urban
    urban_traveltime = p.urbanlength/urbanspeed
    highway_traveltime = p.highwaylength/highwayspeed
    uassign = nothing
    if p.urbanlength <= p.highwaylength
        uassign = highway_traveltime <= urban_traveltime ? p.shortest : 0
    else
        uassign = highway_traveltime <= urban_traveltime ? 1 : (1-p.shortest)
    end

    return uassign
end
immutable DTAInitDist end
sampletype(::Type{DTAInitDist}) = DTAstate
function rand(rng::AbstractRNG, d::DTAInitDist)
    return DTAstate(0, 0, 0.6)
end
initial_state_distribution(::DTAProblem) = DTAInitDist()
#=
function initial_state(p::DTAProblem, rng::AbstractRNG)
    return s0
end
=#
function observation(pp::DTAProblem, s::DTAstate, a::DTAaction, sp::DTAstate)
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    new_highways = sp.highway
    new_urbans = sp.urban
    uassigns = s.uassign
    assigns = a.assign
    new_demands = max(rand(Normal(p.lambda, p.lambda*0.1)), 0)

    qhighway = (p.connected*new_demands.*assigns + (1-p.connected)*new_demands.*uassigns)/p.tau
    k = qhighway/p.highwayspeed
    #if qhighway/p.highwayspeed >= 0.8*p.capacities/p.queuespeed/p.tau
    #    k = 0.8*p.capacities/p.queuespeed/p.tau
    #end
    u = (qhighway - p.capacities/p.tau)/(p.capacities/p.tau/p.queuespeed - k)

    if new_highways == 0
        #assert(u <= 0)
        if u >= 0
            u = -0.1
        end
        tt = s.highway/abs(u)
        highwaymean = highways/2*(tt/p.tau)
    else
        highwaymean = (highways+new_highways)/2
    end
    highwaymeanspeed = p.highwaylength/(highwaymean/p.queuespeed + (p.highwaylength-highwaymean)/p.highwayspeed)
    if a.assign >= 0.4
        hn = Normal(highwaymeanspeed, 1.0)
    else
        hn = Normal(a.assign*highwaymeanspeed+(1-a.assign)*p.highwayspeed, 5.0)
    end

    urban_meandensity = (urbans + new_urbans)/2/p.urbanlength

    urbanmeanspeed = p.maxurbanspeed*(1 - urban_meandensity/p.maxurbandensity)
    if urbanmeanspeed <= 2.5
        urbanmeanspeed = 2.5
    end
    if (1-a.assign) >= 0.4
        un = Normal(urbanmeanspeed, 1.0)
    else
        un = Normal(((1-a.assign)*urbanmeanspeed)+(a.assign*p.maxurbanspeed), 5.0)
    end


    return Records(hn ,un)
end
struct Records
    hn::Normal{Float64}
    un::Normal{Float64}
end

function rand(rng::AbstractRNG, r::Records)
    o_highway = max(rand(rng, r.hn), 0)
    o_urban = max(rand(rng, r.un), 0)
    return DTAobservation(o_highway, o_urban)
end

function pdf(r::Records, o::DTAobservation)
    p = pdf(r.hn, o.highway)*pdf(r.un, o.urban)
    return p
end

function observation(p::DTAProblem, s::DTAstate, a::DTAaction, new_highways::Float64, new_urbans::Float64, u::Float64, rng::AbstractRNG)
    #highwaymean = nothing
    if new_highways == 0
        assert(u <= 0)
        tt = s.highway/abs(u)
        highwaymean = s.highway/2*(tt/p.tau)
    else
        highwaymean = (s.highway+new_highways)/2
    end
    highwaymeanspeed = p.highwaylength/(highwaymean/p.queuespeed + (p.highwaylength-highwaymean)/p.highwayspeed)
    if a.assign >= 0.4
        hn = Normal(highwaymeanspeed, 1.0)
    else
        hn = Normal(a.assign*highwaymeanspeed+(1-a.assign)*p.highwayspeed, 5.0)
    end
    o_highway = max(rand(rng, hn), 0)

    urban_meandensity = (s.urban + new_urbans)/2/p.urbanlength

    urbanmeanspeed = p.maxurbanspeed*(1 - urban_meandensity/p.maxurbandensity)
    if urbanmeanspeed <= 2.5
        urbanmeanspeed = 2.5
    end
    if (1-a.assign) >= 0.4
        un = Normal(urbanmeanspeed, 1.0)
    else
        un = Normal(((1-a.assign)*urbanmeanspeed)+(a.assign*p.maxurbanspeed), 5.0)
    end
    o_urban = max(rand(rng, un), 0)

    return DTAobservation(o_highway, o_urban)
end

function reward(p::DTAProblem, s::DTAstate, a::DTAaction, sp::DTAstate, u::Float64, new_demands::Float64)
    if isterminal(p, sp)
        return -10000.0
    end
    if sp.highway == 0
        assert(u <= 0)
        tt = s.highway/abs(u)
        highwaymean = s.highway/2*(tt/p.tau)
    else
        highwaymean = (s.highway+sp.highway)/2
    end
    cost1 = (highwaymean*p.cost_queue + (p.highwaylength - highwaymean)*p.cost_highway)

    urban_meandensity = (s.urban + sp.urban)/2/p.urbanlength

    urbanmeanspeed = p.maxurbanspeed*(1 - urban_meandensity/p.maxurbandensity)
    if urbanmeanspeed <= 2.5
        urbanmeanspeed = 2.5
    end

    cost2 = p.itpcost[urbanmeanspeed]*p.urbanlength
    #println(cost1,", ", cost2)
    return -(cost1*p.connected*new_demands*a.assign + cost2*p.connected*new_demands*(1-a.assign))*0.000621371
end

discount(p::DTAProblem) = mdp(p).discount
isterminal(p::DTAProblem, s::DTAstate) = mdp(p).highwaylength <= s.highway
immutable ActionSpace end
rand(rng::AbstractRNG, ::ActionSpace) = DTAaction(rand(rng))
actions(::DTAPOMDP) = ActionSpace()

function baseline_action(pp::DTAProblem, s::DTAstate)
    p = mdp(pp)
    cost1 = s.highway*p.cost_queue + (p.highwaylength - s.highway)*p.cost_highway
    urban_meandensity = s.urban/p.urbanlength
    urbanmeanspeed = p.maxurbanspeed*(1 - urban_meandensity/p.maxurbandensity)
    cost2 = p.itpcost[urbanmeanspeed]*p.urbanlength
    if cost1 <= cost2
        return DTAaction(1.0)
    else
        return DTAaction(0.0)
    end
end
###############
@with_kw immutable SymmetricNormalResampler
    n::Int
    std::Float64
end

function ParticleFilters.resample(r::SymmetricNormalResampler, b::WeightedParticleBelief, rng::AbstractRNG)
    collection = resample(LowVarianceResampler(r.n), b, rng)
    ps = particles(collection)
    for i in 1:r.n
        ps[i] += r.std*randn(rng, 2)
    end
    return collection
end

@with_kw immutable MinPopResampler
    n::Int
    min_pop::Int
    std::Float64
end

function ParticleFilters.resample(r::MinPopResampler, b, rng::AbstractRNG)
    collection = resample(LowVarianceResampler(r.n), b, rng)
    ps = particles(collection)
    nu = length(unique(ps))
    if r.min_pop > nu
        is = rand(rng, 1:r.n, r.min_pop - nu)
        for i in is
            ps[i] += r.std*randn(rng, 2)
        end
    end
    @show length(unique(ps))
    @show length(ps)
    return collection
end


immutable ObsAdaptiveParticleFilter{P<:POMDP,S,R,RNG<:AbstractRNG} <: Updater
    pomdp::P
    resample::R
    max_frac_replaced::Float64
    rng::RNG
    _pm::Vector{S}
    _wm::Vector{Float64}
end

function ObsAdaptiveParticleFilter(p::POMDP, resample, max_frac_replaced, rng::AbstractRNG)
    S = state_type(p)
    return ObsAdaptiveParticleFilter(p, resample, max_frac_replaced, rng, S[], Float64[])
end

POMDPs.initialize_belief{S}(up::ObsAdaptiveParticleFilter{S}, d::Any) = resample(up.resample, d, up.rng)
POMDPs.update(up::ObsAdaptiveParticleFilter, b, a, o) = update(up, resample(up.resample, b, up.rng), a, o)

function POMDPs.update(up::ObsAdaptiveParticleFilter, b::ParticleFilters.ParticleCollection, a, o)
    if n_particles(b) > 2*up.resample.n
        b = resample(up.resample, b, up.rng)
    end

    ps = particles(b)
    pm = up._pm
    wm = up._wm
    resize!(pm, 0)
    resize!(wm, 0)

    all_terminal = true
    for i in 1:n_particles(b)
        s = ps[i]
        if !isterminal(up.pomdp, s)
            all_terminal = false
            sp = generate_s(up.pomdp, s, a, up.rng)
            push!(pm, sp)
            od = observation(up.pomdp, s, a, sp)
            push!(wm, pdf(od, o))
        end
    end
    ws = sum(wm)
    if all_terminal || ws < eps(1.0/length(wm))
        # warn("All states in particle collection were terminal.")
        return initialize_belief(up, reset_distribution(up.pomdp, b, a, o))
    end

    pc = resample(up.resample, WeightedParticleBelief{state_type(up.pomdp)}(pm, wm, ws, nothing), up.rng)
    ps = particles(pc)

    mpw = max_possible_weight(up.pomdp, a, o)
    frac_replaced = up.max_frac_replaced*max(0.0, 1.0 - maximum(wm)/mpw)
    n_replaced = floor(Int, frac_replaced*length(ps))
    is = randperm(up.rng, length(ps))[1:n_replaced]
    for i in is
        ps[i] = new_particle(up.pomdp, b, a, o, up.rng)
    end
    return pc
end

# POMCPOW.belief_type(::Type{ObsAdaptiveParticleFilter{Vec2}}, ::Type{LightDark2DTarget}) = POWNodeBelief{Vec2, Vec2, Vec2, LightDark2DTarget}



reset_distribution(p::POMDP, b, a, o) = initial_state_distribution(p)




max_possible_weight(pomdp::DTAPOMDP, a::DTAaction, o::DTAobservation) = 0.0

new_particle(pomdp::DTAPOMDP, a::DTAaction, o::DTAobservation) = error("shouldn't get here")



end
