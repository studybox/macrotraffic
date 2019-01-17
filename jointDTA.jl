module jointDTA

using macrotraffic
using POMDPs
using POMDPToolbox
using ParticleFilters

using Distributions
using StaticArrays
using Parameters
using Interpolations

using PyCall
include("/media/mSATA/UM/Simulation/constants.jl")
include("/media/mSATA/UM/Upper routing simulation/utils.jl")
@pyimport sys
@pyimport os
unshift!(PyVector(pyimport("sys")["path"]), os.path[:join](SUMO_HOME, "tools"))
unshift!(PyVector(pyimport("sys")["path"]), "/media/mSATA/UM/Simulation")
@pyimport sumolib
M2MILE = 0.000621371
ANNARBOR = sumolib.net[:readNet]("/media/mSATA/UM/Simulation/f_AnnArbor1.6.net.xml", withPrograms=true)
Zones = read_zonefile("/media/mSATA/UM/Simulation/finalclusters.txt")
sumomap = load_sumo_map2(Zones[1], Zones[2]);
linkname = Dict()
for (idx, ln) in enumerate(data_links[:link])
    linkname[ln] = idx
end
M2MILE = 0.000621371
CARLENGTH = 7.5
R1 = ["10272-1", "10272-1-Pocket","10284-1","10284-1-Pocket","10352-1","10352-1-Pocket","10359-1","10359-1-Pocket",
      "10360-1","10387-1","10462-1","10462-1-Pocket","10444-1","10444-1-Pocket","10577-1","10491-1","10491-1-Pocket",
      "10631-1","10769-1","10789-1","10789-1-Pocket","10709-1","13656-1","13657-1","11027-1","11035-1","11035-1-Pocket",
      "11082-1","11082-1-Pocket","11247-1","11247-1-Pocket"]

R2 = ["10324-0","12491-0","12499-1","12551-1","12557-1","12602-1","12602-1-Pocket","12537-1",
      "12543-1","13600-1","13597-1","13597-1-Pocket","12680-0","11336-1","11334-1","11328-1",
      "11328-1-Pocket","11215-1","11215-1-Pocket","11187-1","11187-1-Pocket","11254-0","11254-0-Pocket"]
r = vcat(R1, R2)
nns = Set()
for l in r
    tonodeid = NETEDGESLOOKUP[l][:getToNode]()[:getID]()
    fromnodeid = NETEDGESLOOKUP[l][:getFromNode]()[:getID]()
    push!(nns, fromnodeid)
    push!(nns, tonodeid)
end
n = collect(nns)

t = []


for i = 1:length(r)
    ln = NETEDGES[r[i]][2]
    lid = linkname[ln]
    typee = data_links[:type][lid]
    push!(t, typee)
end
information = Dict("edges"=>r, "nodes"=>n, "etypes"=>t)
NET = import_network(ANNARBOR, information)
add_route!(NET, R1)
add_route!(NET, R2)

DEFAULT_URBAN_TIME = compute_default_time(NET, NET.routes[1].id)
DEFAULT_HIGHWAY_TIME = compute_default_time(NET, NET.routes[2].id)
import ParticleFilters: obs_weight

importall POMDPs
import Base: rand

export DTAstate, DTAaction, DTAobservation, DTAMDP, DTAPOMDP, mdp, obs_weight, ObsAdaptiveParticleFilter, generate_sr, isterminal, baseline_action, baseline_action2, actions, generate_actions, generate_orders

struct DTAstate
    step::Int
    highway::Float64
    urban::Float64
end

struct DTAjoint
    e::Int
    nodeidx::Int
    #assign :: Vector{Float64}
end
const DTAaction = Vector{DTAjoint}
#=
struct DTAaction
    assign :: Float64
end
=#
struct DTAobservation
    highway :: Float64
    urban :: Float64
end

immutable DTAMDP <: MDP{DTAstate, DTAaction}
    maxstep::Int # = 13
    actiondev::Int
    lambdas::Vector{Float64}# = [0.8]
    fastest::Float64# = 0.5
    connected::Float64# = 0.5
    tau::Int# = 5*60
    omega::Float64# = 40
    discount::Float64# = 0.999999999
    E::Int# = 1
    numagents::Int# = 1
    joint_actions::Vector{DTAaction}
    joint_action_groups::Vector{Vector{Int}}
    joint_state_groups::Vector{Dict{Vector{Int},Vector{Int}}}
    elimination_order::Vector{Int}

    function DTAMDP(;maxstep::Int=13, actiondev::Int=11,
                    lambdas::Vector{Float64}=[0.8,0.8,0.8],
                    fastest::Float64=0.5,
                    connected::Float64=0.5,
                    tau::Int=300,
                    omega::Float64=40.0,
                    discount::Float64=0.9999999,
                    E::Int=2, numagents::Int=3,
                    joint_action_groups::Vector{Vector{Int}}=[[1,2],[2,3]],
                    joint_state_groups::Vector{Dict{Vector{Int},Vector{Int}}}=[Dict([1,1]=>[1,2,3],[2,2]=>[3,4,5]),Dict([1,1]=>[2,3,4],[2,2]=>[4,5,6])],
                    elimination_order::Vector{Int}=[1,3,2])
        joint_actions = generate_joint_actions(joint_action_groups, E, numagents, actiondev)
        return new(
                maxstep,
                actiondev,
                lambdas,
                fastest,
                connected,
                tau,
                omega,
                discount,
                E, numagents,
                joint_actions,
                joint_action_groups,
                joint_state_groups,
                elimination_order
        )
    end
end

@with_kw immutable DTAPOMDP <: POMDP{DTAstate, DTAaction, DTAobservation}
    mdp ::DTAMDP = DTAMDP()
    obs_std::Float64 = 0.5
end

const DTAProblem = Union{DTAMDP, DTAPOMDP}



mdp(p::DTAPOMDP) = p.mdp
mdp(p::DTAMDP) = p

function generate_joint_actions(joint_action_groups, E, numagents, actiondev)
    joint_actions = []
    for i = 1:E
        jj = []
        na  = length(joint_action_groups[i])
        for j = 1:actiondev^na
            ja = DTAjoint(i, j)
            push!(jj, ja)
        end
        push!(joint_actions, jj)
    end
    return joint_actions
end
function generate_actions(pp::DTAProblem)
    p = mdp(pp)
    A = DTAaction[]
    dimensions = Tuple([0:p.actiondev-1 for i=1:p.numagents])
    for ii = 1:p.actiondev^p.numagents
        a = ind2sub(dimensions, ii)
        #println(a)
        a = [i for i in a]
        jas = DTAjoint[]
        for (e,j) in enumerate(p.joint_action_groups)
            #println(j)
            ja = a[j]
            #println(ja)
            d = Tuple([0:p.actiondev-1 for i=1:length(j)])
            #println(d)
            nodeidx = sub2ind(d,ja...)
            #println(nodeidx)
            joint = DTAjoint(e,nodeidx)
            push!(jas, joint)
        end
        push!(A, jas)
    end
    A
end

function generate_orders(pp::DTAProblem)
    p = mdp(pp)
    eo = p.elimination_order
    groups = deepcopy(p.joint_action_groups)
    used = zeros(length(groups))
    retO = []
    retJ = []
    for (idxo, o) in enumerate(eo)
        newjointid = length(groups)+1
        agentids = []
        jointids = []
        for (idxaa, aa) in enumerate(groups)
            if used[idxaa] == 0 && o in aa
                append!(agentids, aa)
                push!(jointids, idxaa)
                used[idxaa] = 1
            end
        end
        agentids = sort(unique(agentids))
        #println(agentids)
        newjointgroup = filter(e->e≠o, agentids)
        #println(agentids)
        push!(groups, newjointgroup)
        push!(used, 0)
        dim = Tuple([0:p.actiondev-1 for i=1:length(agentids)])
        #println(dim)
        Odict = Dict()
        Jlist = []
        for ii = 1:p.actiondev^length(agentids)
            if length(dim) > 1

                agentacts = ind2sub(dim, ii)
                agentacts = [i for i in agentacts]

            else
                agentacts = [ii-1]
            end
            jointacts = []
            nodeidxs = []
            for e in jointids
                jindexs = []
                for ag in groups[e]
                    push!(jindexs, findfirst(agentids, ag))
                end
                ja = agentacts[jindexs]
                d = Tuple([0:p.actiondev-1 for i=1:length(jindexs)])
                if length(jindexs) > 1
                    nodeidx = sub2ind(d,ja...)
                else
                    nodeidx = ja[1]+1
                end
                push!(nodeidxs, nodeidx)
                joint = DTAjoint(e,nodeidx)
                push!(jointacts, joint)
            end
            push!(Jlist, jointacts)
            newjindexs = []
            for ag in newjointgroup
                push!(newjindexs, findfirst(agentids, ag))
            end
            newja = agentacts[newjindexs]
            if length(newjindexs) > 1
                newd = Tuple([0:p.actiondev-1 for i=1:length(newjindexs)])
                newnodeidx = sub2ind(newd, newja...)

            elseif length(newjindexs) == 1
                newnodeidx = newja[1]+1

            else
                newnodeidx = -jointids[1]
            end

            #newjoint = DTAjoint(newjointid, newnodeidx)
            #Odict[jointacts] = (newjointid, newnodeidx)
            Odict[nodeidxs] = (newjointid, newnodeidx)
        end
        push!(retO, Odict)
        push!(retJ, Jlist)
    end
    return (retO, retJ, groups)
end





LHIGH1_1 = 278.71
LHIGH1_2 = 514.99
LHIGH2 = 5324.61
LHIGH3 = 402.35
LHIGH4 = 1048.9

CHIGH1_1 = 83.613
CHIGH1_2 = 51.50
CHIGH2 = 670.22
CHIGH3 = 40.235

LURBAN1 = 1000
LURBAN2 = 5389.79

CURBAN1 = 220
CURBAN2 = 1156.509

VHIGH = 2.73308349191062*3
QHIGH = 0.19715888888888888*3

VURBAN = 1.057172133849443
QURBAN = 0.20456333333333332

function get_uassigns(p, highway, urban)
    highway_time = 0.0
    highway_time += LHIGH4/20.12
    if highway >= CHIGH3
        highway_time += LHIGH3/VHIGH
        resh = highway - CHIGH3
        if resh >= CHIGH2
            highway_time += LHIGH2/VHIGH
            resh = resh - CHIGH2
            if resh >= CHIGH1_2
                highway_time += LHIGH1_2/VHIGH
                resh = resh - CHIGH1_2
                highway_time += resh*CARLENGTH/(3*VHIGH)
                highway_time += max(LHIGH1_1 - resh*CARLENGTH/3, 0)/20.12
            else
                highway_time += resh*CARLENGTH/VHIGH
                highway_time += max(LHIGH1_2 - resh*CARLENGTH, 0)/31.29
                highway_time += LHIGH1_1/20.12
            end
        else
            highway_time += resh*CARLENGTH/(1.25*VHIGH)
            highway_time += max(LHIGH2 - resh*CARLENGTH/1.25, 0)/31.29
            highway_time += LHIGH1_2/31.29
            highway_time += LHIGH1_1/20.12
        end
    else
        highway_time += highway*CARLENGTH/VHIGH
        highway_time += max(LHIGH3 - highway*CARLENGTH, 0)/17.88
        highway_time += LHIGH2/31.29
        highway_time += LHIGH1_2/31.29
        highway_time += LHIGH1_1/20.12
    end
###############

    urban_time = 0.0
    urban_time += LURBAN2/17.88
    urban_time += urban*CARLENGTH/(2.2*VURBAN)
    urban_time += max(LURBAN1 - (urban*CARLENGTH/2.2), 0)/20.12

    if urban_time < highway_time
        return 0.0
    else
        return p.fastest
    end
end

function compute_time(s::DTAstate)
    highway = s.highway
    urban = s.urban

    highway_time = 0.0
    highway_time += LHIGH4/20.12
    if highway >= CHIGH3
        highway_time += LHIGH3/VHIGH
        resh = highway - CHIGH3
        if resh >= CHIGH2
            highway_time += LHIGH2/VHIGH
            resh = resh - CHIGH2
            if resh >= CHIGH1_2
                highway_time += LHIGH1_2/VHIGH
                resh = resh - CHIGH1_2
                highway_time += resh*CARLENGTH/(3*VHIGH)
                highway_time += max(LHIGH1_1 - resh*CARLENGTH/3, 0)/20.12
            else
                highway_time += resh*CARLENGTH/VHIGH
                highway_time += max(LHIGH1_2 - resh*CARLENGTH, 0)/31.29
                highway_time += LHIGH1_1/20.12
            end
        else
            highway_time += resh*CARLENGTH/(1.25*VHIGH)
            highway_time += max(LHIGH2 - resh*CARLENGTH/1.25, 0)/31.29
            highway_time += LHIGH1_2/31.29
            highway_time += LHIGH1_1/20.12
        end
    else
        highway_time += highway*CARLENGTH/VHIGH
        highway_time += max(LHIGH3 - highway*CARLENGTH, 0)/17.88
        highway_time += LHIGH2/31.29
        highway_time += LHIGH1_2/31.29
        highway_time += LHIGH1_1/20.12
    end
###############

    urban_time = 0.0
    urban_time += LURBAN2/17.88
    urban_time += urban*CARLENGTH/(2.2*VURBAN)
    urban_time += max(LURBAN1 - (urban*CARLENGTH/2.2), 0)/20.12

    return (highway_time, urban_time)
end

function get_new_demands(lambda, T)
    u = Poisson(T*lambda)
    newveh = rand(u)
    return newveh
end




function generate_s(pp::DTAPOMDP, s::DTAstate, a::DTAaction, rng::AbstractRNG)
    #println("start s")
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    assigns = a.assign #min(max(a.assign, 0.0), 1.0)
    uassigns = get_uassigns(p, highways, urbans)
    new_demands = get_new_demands(p.lambda, p.tau)
    Qhighway = (p.connected*new_demands.*assigns) + ((1-p.connected)*new_demands.*uassigns)
    Qurban = (p.connected*new_demands.*(1-assigns)) + ((1-p.connected)*new_demands.*(1-uassigns))

    new_highways = max(highways - QHIGH*p.tau + Qhighway, 0)
    new_urbans = max(urbans - QURBAN*p.tau + Qurban, 0)

    sp = DTAstate(s.step+1, new_highways, new_urbans)
    #dd = generate_demands([(NET.routes[1].O,NET.routes[1].D), (NET.routes[2].O,NET.routes[2].D)], [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns], NET.routes)
    #queuedict = Dict()
    #queuedict[8] = s.urban
    #queuedict[44] = s.highway
    #urate = Dict()
    #urate[p.lambda*(1.0-p.connected)] = (1, 2, p.shortest)


    #newnet, cost = state_step!(NET, dd, queuedict, urate, rng, [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns])
    #=
    urban_queue = 0.0
    urban_queue += dd[(newnet.routes[1].O,newnet.routes[1].roads[1])][1].vehicles
    roadIds = newnet.routes[1].roads
    for (idx, rid) in enumerate(roadIds)
        road = get_road(newnet, rid)
        if idx == length(roadIds)
            nexto = -1
        else
            nexto = roadIds[idx+1]
        end
        flow = nothing
        for f in road.collection[nexto]
            if f.routeid == newnet.routes[1].id
                flow = f
                break
            end
        end
        urban_queue += flow.queue
        if rid == 8
            break
        end
    end

    highway_queue = 0.0
    highway_queue += dd[(newnet.routes[2].O,newnet.routes[2].roads[1])][1].vehicles
    roadIds = newnet.routes[2].roads
    for (idx, rid) in enumerate(roadIds)
        road = get_road(newnet, rid)
        if idx == length(roadIds)
            nexto = -1
        else
            nexto = roadIds[idx+1]
        end
        flow = nothing
        for f in road.collection[nexto]
            if f.routeid == newnet.routes[2].id
                flow = f
                break
            end
        end
        highway_queue += flow.queue
        if rid == 44
            break
        end
    end
    sp = DTAstate(s.step+1, highway_queue, urban_queue)
    =#
    return sp
end

function compute_energy(s::DTAstate)
    highway = s.highway
    highway_energy = 0.0
    highway_energy += (LHIGH4-0)*compute_energy_rate(20.12, 20.12)*M2MILE
    if highway >= CHIGH3
        highway_energy += LHIGH3*compute_energy_rate(VHIGH, 17.88)*M2MILE
        resh = highway - CHIGH3
        if resh >= CHIGH2
            highway_energy += LHIGH2*compute_energy_rate(VHIGH, 31.29)*M2MILE
            resh = resh - CHIGH2
            if resh >= CHIGH1_2
                highway_energy += LHIGH1_2*compute_energy_rate(VHIGH, 17.88)*M2MILE
                resh = resh - CHIGH1_2
                highway_energy += resh*CARLENGTH*compute_energy_rate(VHIGH, 20.12)*M2MILE/3
                highway_energy += max(LHIGH1_1 - resh*CARLENGTH/3, 0)*compute_energy_rate(20.12, 20.12)*M2MILE
            else
                highway_energy += resh*CARLENGTH*compute_energy_rate(VHIGH, 17.88)*M2MILE
                highway_energy += max(LHIGH1_2 - resh*CARLENGTH)*compute_energy_rate(31.29, 17.88)*M2MILE
                highway_energy += LHIGH1_1*compute_energy_rate(20.12, 20.12)*M2MILE
            end
        else
            highway_energy += resh*CARLENGTH*compute_energy_rate(VHIGH, 31.29)*M2MILE/1.25
            highway_energy += max(LHIGH2 - resh*CARLENGTH/1.25, 0)*compute_energy_rate(31.29, 31.29)*M2MILE
            highway_energy += LHIGH1_2*compute_energy_rate(31.29, 17.88)*M2MILE
            highway_energy += LHIGH1_1*compute_energy_rate(20.12, 20.12)*M2MILE
        end
    else
        highway_energy += highway*CARLENGTH*compute_energy_rate(VHIGH, 17.88)*M2MILE
        highway_energy += max(LHIGH3 - highway*CARLENGTH, 0)*compute_energy_rate(17.88, 17.88)*M2MILE
        highway_energy += LHIGH2*compute_energy_rate(31.29, 31.29)*M2MILE
        highway_energy += LHIGH1_2*compute_energy_rate(31.29, 17.88)*M2MILE
        highway_energy += LHIGH1_1*compute_energy_rate(20.12, 20.12)*M2MILE
    end

    urban = s.urban
    urban_energy = 0.0
    urban_energy += LURBAN2*compute_energy_rate(17.88, 17.88)*M2MILE
    urban_energy += (urban*CARLENGTH/2.2)*compute_energy_rate(VURBAN, 20.12)*M2MILE
    urban_energy += max(LURBAN1 - (urban*CARLENGTH/2.2), 0)*compute_energy_rate(VURBAN, 20.12)*M2MILE
    return (highway_energy, urban_energy)
end
function reward(pp::DTAPOMDP, s::DTAstate, a::DTAaction, sp::DTAstate)
    p = mdp(pp)
    #=
    if sp.step >= p.maxstep
        return 0
    end
    =#
    assigns = a.assign #min(max(a.assign, 0.0), 1.0)
    highways = s.highway
    urbans = s.urban
    uassigns = get_uassigns(p, highways, urbans)
    rhighway = (p.connected*assigns) + ((1-p.connected)*uassigns)
    rurban = (p.connected*(1-assigns)) + ((1-p.connected)*(1-uassigns))
    energy_s = compute_energy(s)
    energy_sp = compute_energy(sp)

    cost = ((energy_s[1]+energy_sp[1])*0.5*rhighway) + ((energy_s[2]+energy_sp[2])*0.5*rurban)

    return 10.0 - cost


end
function generate_sr(pp::DTAPOMDP, s::DTAstate, a::DTAaction, rng::AbstractRNG)
    #println("start sr")
    p = mdp(pp)
    highways = s.highway
    urbans = s.urban
    assigns = a.assign #min(max(a.assign, 0.0), 1.0)
    uassigns = get_uassigns(p, highways, urbans)
    new_demands = get_new_demands(p.lambda, p.tau)
    Qhighway = (p.connected*new_demands.*assigns + (1-p.connected)*new_demands.*uassigns)
    Qurban = (p.connected*new_demands.*(1-assigns) + (1-p.connected)*new_demands.*(1-uassigns))

    new_highways = max(highways - QHIGH*p.tau + Qhighway, 0)
    new_urbans = max(urbans - QURBAN*p.tau + Qurban, 0)

    sp = DTAstate(s.step+1, new_highways, new_urbans)

    r = reward(pp, s, a, sp)
    #=
    dd = generate_demands([(NET.routes[1].O,NET.routes[1].D), (NET.routes[2].O,NET.routes[2].D)], [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns], NET.routes)
    queuedict = Dict()
    queuedict[8] = s.urban
    queuedict[44] = s.highway
    urate = Dict()
    urate[p.lambda*(1.0-p.connected)] = (1, 2, p.shortest)
    newnet, cost = state_step!(NET, dd, queuedict, urate, rng, [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns])

    urban_queue = 0.0
    urban_queue += dd[(newnet.routes[1].O,newnet.routes[1].roads[1])][1].vehicles
    roadIds = newnet.routes[1].roads
    for (idx, rid) in enumerate(roadIds)
        road = get_road(newnet, rid)
        if idx == length(roadIds)
            nexto = -1
        else
            nexto = roadIds[idx+1]
        end
        flow = nothing
        for f in road.collection[nexto]
            if f.routeid == newnet.routes[1].id
                flow = f
                break
            end
        end
        urban_queue += flow.queue
        if rid == 8
            break
        end
    end

    highway_queue = 0.0
    highway_queue += dd[(newnet.routes[2].O,newnet.routes[2].roads[1])][1].vehicles
    roadIds = newnet.routes[2].roads
    for (idx, rid) in enumerate(roadIds)
        road = get_road(newnet, rid)
        if idx == length(roadIds)
            nexto = -1
        else
            nexto = roadIds[idx+1]
        end
        flow = nothing
        for f in road.collection[nexto]
            if f.routeid == newnet.routes[2].id
                flow = f
                break
            end
        end
        highway_queue += flow.queue
        if rid == 44
            break
        end
    end
    sp = DTAstate(s.step+1, highway_queue, urban_queue)
    #r = reward(p, s, a, sp, u, new_demands)
    r = 1000-cost
    #println("end sr")
    #println(actions(pp))
    #println(assigns)

    #if assigns == 0.0
#        println(r)
    #end
    #if assigns == 0.9 || assigns == 1.0
    #    r = -10000.0
    #end
    =#

    return (sp, r)
end


immutable DTAInitDist end
sampletype(::Type{DTAInitDist}) = DTAstate
function rand(rng::AbstractRNG, d::DTAInitDist)
    return DTAstate(0, 0.0, 0.0)
end
initial_state_distribution(::DTAProblem) = DTAInitDist()
#=
function initial_state(p::DTAProblem, rng::AbstractRNG)
    return s0
end
=#
function observation(pp::DTAProblem, s::DTAstate, a::DTAaction, sp::DTAstate)
    #println("start observe")
    p = mdp(pp)
    rng = nothing

    assigns = a.assign #min(max(a.assign, 0.0), 1.0)
    #dd = generate_demands([(NET.routes[1].O,NET.routes[1].D), (NET.routes[2].O,NET.routes[2].D)], [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns], NET.routes)
    #queuedict = Dict()
    #queuedict[8] = sp.urban
    #queuedict[44] = sp.highway

    #time_urban = sp.urban #compute_time2(NET, queuedict, NET.routes[1].id)
    #time_highway = sp.highway #compute_time2(NET, queuedict, NET.routes[2].id)
    #dtime_urban = DEFAULT_URBAN_TIME
    #dtime_highway = DEFAULT_HIGHWAY_TIME
    hn = Normal(s.highway, 0.0001)
    un = Normal(s.urban, 0.0001)
    #=
    if a.assign >= 0.4
        hn = Normal(time_highway, 10.0)
    else
        hn = Normal(a.assign*time_highway+(1-a.assign)*dtime_highway, 60.0)
    end

    if (1.0-a.assign) >= 0.4
        un = Normal(time_urban, 10.0)
    else
        un = Normal(((1.0-a.assign)*time_urban)+(a.assign*dtime_urban), 60.0)
    end
    =#
    #println("end observe")
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



discount(p::DTAProblem) = mdp(p).discount
isterminal(p::DTAProblem, s::DTAstate) = s.highway >= 500# mdp(p).maxstep
immutable ActionSpace end

#rand(rng::AbstractRNG, ::ActionSpace) = DTAaction(rand(rng))
#actions(::DTAPOMDP) = ActionSpace()
n_actions(::DTAPOMDP) = 11
actions(::DTAPOMDP) = [DTAaction(x) for x = 0.0:0.1:1.0]

function baseline_action(pp::DTAProblem, s::DTAstate)
    p = mdp(pp)
    rng = nothing
    #NET = import_network(ANNARBOR, information)
    #add_route!(NET, R1)
    #add_route!(NET, R2)
    #=
    assigns = 0.5
    dd = generate_demands([(NET.routes[1].O,NET.routes[1].D), (NET.routes[2].O,NET.routes[2].D)], [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns], NET.routes)
    queuedict = Dict()
    queuedict[8] = s.urban
    queuedict[44] = s.highway
    urate = Dict()
    urate[p.lambda*(1.0-p.connected)] = (1, 2, p.shortest)
    state =  create_state!(NET, dd, queuedict, urate, rng)
    #energy_urban = compute_energy_cost2(NET, queuedict, NET.routes[1].id)
    #energy_highway = compute_energy_cost2(NET, queuedict, NET.routes[2].id)
    energy_urban = compute_energy_cost(state, dd, NET.routes[1].id)
    energy_highway = compute_energy_cost(state, dd, NET.routes[2].id)
    =#
    energy_highway, energy_urban = compute_energy(s)

    if energy_urban <= energy_highway #&& s.urban <= 750
        return DTAaction(0.0)
    else
        return DTAaction(1.0)
    end
end
function baseline_action2(pp::DTAProblem, s::DTAstate)
    p = mdp(pp)
    rng = nothing
    #NET = import_network(ANNARBOR, information)
    #add_route!(NET, R1)
    #add_route!(NET, R2)
    #=
    assigns = 0.5
    dd = generate_demands([(NET.routes[1].O,NET.routes[1].D), (NET.routes[2].O,NET.routes[2].D)], [p.lambda*p.connected*(1.0-assigns), p.lambda*p.connected*assigns], NET.routes)
    queuedict = Dict()
    queuedict[8] = s.urban
    queuedict[44] = s.highway
    urate = Dict()
    urate[p.lambda*(1.0-p.connected)] = (1, 2, p.shortest)
    state =  create_state!(NET, dd, queuedict, urate, rng)
    =#
    time_highway, time_urban = compute_time(s)
    #time_urban = compute_time(state, dd, NET.routes[1].id)
    #time_highway = compute_time(state, dd, NET.routes[2].id)

    if time_urban <= time_highway #&& s.urban <= 750
        return DTAaction(0.0)
    else
        return DTAaction(1.0)
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
