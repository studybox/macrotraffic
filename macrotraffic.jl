module macrotraffic

using PyCall
using DataStructures
using Distributions
using StaticArrays
using Parameters
using Graphs
import Base: show
export Flow, Urban, Junction, Connection, Demand, Network, import_network
CARLENGTH = 10
Tu = 2 #sec
@with_kw mutable struct Flow
    routeid::Int
    roadindex::Int = 0
    queue::Int = 0
    arrival1::Int = 0
    arrival2::DefaultDict = DefaultDict{Int, Int}(0)
    depart::Int = 0
end

mutable struct Urban
    id::Int
    storage::Int
    collection::Dict{Int, Array{Flow,1}}
    maxstorage::Int
    delay::Int
    spdlim::Float64
    length::Float64
    numlane::Int
end
struct Route
    id::Int
    roads::Array{Int, 1}
    O::Int
    D::Int
end
#=
struct Freeway
end
struct Ramp
end
=#
struct Connection
    fromEdge::Int
    toEdge::Int
    throughlane::Int
    satflow::Float64
    direction::String
    red::Float64
end

struct Junction
    id::Int
    name::String
    jtype::String
    inroads::Array{Int, 1}
    outroads::Array{Int, 1}
    connections::Dict{Tuple{Int, Int}, Connection}
    cycle::Float64
end

@with_kw mutable struct Demand
    rate::Int
    O::Int
    D::Int
    routeid::Int
    vehicles::Int = 0
end




const Road = Union{Urban}#, Freeway, Ramp}

function get_source(demands::DefaultDict{Tuple{Int, Int}, Array{Demand,1}}, junctionid::Int, outroadid::Int)
    return demands[(Junctionid, outroadid)]
end

struct Network
    roads::SVector
    junctions::SVector
    routes::Array{Route, 1}
end
Network(rs::Array, js::Array) = Network(SVector{length(rs)}(rs), SVector{length(js)}(js), Array{Int, 1}[])
function get_road(net::Network, id::Int)
    return net.roads[id]
end
function get_flow1(road, did, Qcap)
    flows = road.collection[did]
    q = 0
    for f in flows
        q += f.queue + f.arrival[0]
    end
    return minimum(q, Tu*Qcap)
end
function get_capacity(junction::Junction, oid::Int, did::Int)
    con = junction.connections[(oid, did)]
    p = con.satflow # veh/hr
    r = con.red
    c = junction.cycle
    ra = r+10 >= c ? r : r+10
    Q = p*(c-ra)/c
    return Q
end
function create_vehicles!(demands::Array{Demand, 1}, rng::AbstractRNG)
    pq = Priorityqueue{Int, Int}()
    totalf = 0
    for (idx, d) in enumerate(demands)
        f = create_vehicles!(d, rng)
        enqueue!(pq, idx=>f)
        totalf += f
    end
    return (f, totalf)
end
function create_vehicles!(demand::Demand, rng::AbstractRNG)
    u = Poisson(Tu*demand.rate)
    newveh = rand(u, rng)
    demand.vehicles += newveh
    return demand.vehicles
end

function update_storage!(road::Road, fs::PriorityQueue, tfs::Int)

    if tf <= road.storage
        road.storage -= tf
        return fs
    else
        pq = Priorityqueue{Int, Int}()
        N = length(fs)
        for i = 1:N
            (idx, f) = peek(fs)
            if f <= road.storage/(N-i+1)
                road.storage -= f
                enqueue!(pq, idx=>f)
            else
                q = Int(floor(road.storage/(N-i+1)))
                for (ii, idx2) in enumerate(keys(fs))
                    if ii != N-i+1
                        enqueue!(pq, idx2=>q)
                    else
                        qq = road.storage - q*(N-i)
                        enqueue!(pq, idx2=>qq)
                    end
                end
                road.storage = 0
                break
            end
            dequeue!(fs)
        end
        return pq
    end
end

function update_vehicles!(demand::Array{Demand, 1}, f::PriorityQueue)
    for (idx, d) in enumerate(demand)
        d.vehicles -= f[idx]
        assert(d.vehicles >= 0)
    end
end

function update_arrival1_and_depart!(road, flows, net)
    id = road.id
    for inid in keys(flows)
        f = flows[inid] # the total number of cars in arrival
        inroad = get_road(net, inid)
        ff = inroad.collection[id]

        N = length(ff)
        res = f % N
        share = Int(floor(f / N))
        for (idx, flow) in enumerate(ff)
            arrival1 = idx <= res ? share+1 : share
            inroad.collection[id][idx].depart = arrival1

            routeid = flow.routeid
            roadindex = flow.roadindex
            route = get_route(net, routeid)
            update_arrival1!(road, route, routeid, roadindex, arrival1)
        end



        #update_depart!(inroad, f)
        #update_arrival1!(road, ff, f, net)
    end

end

function update_arrival1!(road, route, routeid, roadindex, arrival1)
    if length(route) == roadindex
        nexto = - 1
    else
        nexto = route[roadindex+1]
    end
    if nexto in road.collection
        hasroute = false

        for flow2 in road.collection[nexto]
            if flow2.routeid == routeid
                flow2.arrival1 = arrival1
                hasroute = true
                break
            end
        end
        if !hasroute
            push!(road.collection[nexto], Flow(routeid = routeid, roadindex = roadindex+1, arrival1=arrival1))
        end
    else
        road.collection[nexto] = [Flow(routeid = routeid, roadindex = roadindex+1, arrival1=arrival1)]
    end
end

function update_arrival!(road, demand, ff, net)
    for (idx, d) in enumerate(demand)
        arrival1 = ff[idx]
        routeid = d.routeid
        roadindex = 1
        route = get_route(net, routeid)
        update_arrival1!(road, route, routeid, roadindex, arrival1)
    end
end

function update_storage!(road)
    for d in keys(road.collection)
        for f in road.collection[d]
            road.storage += f.depart
        end
    end
    road.storage = minimum(road.storage, road.maxstorage)
end

function update_queue!(road)
    for d in keys(road.collection)
        for f in road.collection[d]
            f.queue = max(f.queue + f.arrival2[0] - f.depart , 0)
        end
    end
end
function update_arrival2!(road)
    for d in keys(road.collection)
        for f in road.collection[d]
            f.arrival2[road.delay] = f.arrival2[road.delay] + f.arrival1
            # update to next time step
            f.arrival1 = 0
            new_arrival2 = DefaultDict{Int, Int}(0)
            for t in keys(f.arrival2)
                new_arrival2[t-1] = f.arrival2[t]
            end
            f.arrival2 = new_arrival2
        end
    end
    road.delay = Int(ceil(road.storage*CARLENGTH/(road.spdlim*Tu)))
end
function is_connected(junction::Junction, inid::Int, outid::Int)
    return (inid, outid) in keys(junction.connections)
end
function generate_demands(ODs, rates, routes)
    dd = DefaultDict{Tuple{Int, Int}, Array{Demand, 1}}(Demand[])
    for (idx, od) in enumerate(ODs)
        d = Demand(rate=rates[idx], O=od[1], D=od[2], routeid=routes[idx].id)
        if (od[1], routes[idx].roads[1]) in keys(dd)
            push!(dd[(od[1], routes[idx].roads[1])] , d)
        else
            dd[(od[1], routes[idx].roads[1])] = [d]
        end
    end
    return dd
end
function step(network, demands, rng)
    for junction in network.junctions
        for outid in junction.outroads
            outroad = get_road(network, outid)
            demand = get_source(demands, junction.id, outid)
            if length(demand) != 0
                # compute intentions
                flows1, totalflows1 = create_vehicles!(demand, rng)
                # Must satisfy demand first
                # Compute and update remaining storage
                flows2 = update_storage!(outroad, flows1, totalflows1)
                update_arrival1!(outroad, demand, flows2, net)
                # update demand if not all satisfied
                update_vehicles!(demand, flows2)
            end
            flows1 = Priorityqueue{Int, Int}()
            totalflows1 = 0
            for inid in junction.inroads
                if is_connected(junction.connections, inid, outid)
                    inroad = get_road(network, inid)
                    Qcap = get_capacity(junction, inid, outid)
                    flow1 = get_flow1(inroad, outid, Qcap)
                    enqueue!(flows1, inid=>flow1)
                    totalflows1 += flow1
                end
            end
            # decide actual flow
            flows2 = update_storage!(outroad, flows1, totalflows1)
            update_arrival1_and_depart!(outroad, flows2, network)
        end
    end
    for road in network.roads
        # Update storage
        update_storage!(road)

        update_queue!(road)

        update_arrival2!(road)

    end
end
function add_route!(net::Network, route::Route)
    push!(net.routes, route)
end
function build_network(froms, tos, edgeattr, nodeattr)
    js = []
    for nodeidx = 1:length(nodeattr)
        outroads = nodeattr[nodeidx]["outroad"]
        inroads = nodeattr[nodeidx]["inroad"]
        cons = nodeattr[nodeidx]["connection"]
        ps = nodeattr[nodeidx]["saturation"]
        c = nodeattr[nodeidx]["cycle"]
        rs = nodeattr[nodeidx]["red"]
        j = Junction(nodeidx, inroads, outroads, cons, ps, c, rs)
        push!(js, j)
    end
    es = []
    for edgeidx = 1:length(froms)
        if edgeattr[edgeidx]["type"] == "urban"
            L = edgeattr[edgeidx]["length"]
            numlane = edgeattr[edgeidx]["lane"]
            storage = Int(floor(L * numlane / CARLENGTH))
            junction = js[tos[edgeidx]]
            co = Dict{Int, Array{Flow,1}}()
            for outid in junction.outroads
                if is_connected(junction.connections, edgeidx, outid)
                    co[outid] = Flow[]
                end
            end
            spdlim = edgeattr[edgeidx]["spdlim"]
            e = Urban(edgeidx, storage, co, storage, 1, spdlim, L, numlane)
            push!(es, e)
        end
    end
    return Network(js, es)
end
function import_network(net::PyObject, info)
    edgeids = info["edges"]
    nodeids = info["nodes"]
    etypes = info["etypes"]
    n2id = Dict()
    for idx = 1:length(nodeids)
        n2id[nodeids[idx]] = idx
    end
    l2id = Dict()
    for idx = 1:length(edgeids)
        l2id[edgeids[idx]] = idx
    end
    js = []
    for (index, nodeid) in enumerate(nodeids)
        node = net[:getNode](nodeid)
        incos = node[:getIncoming]()
        inroads = []
        for i in incos
            l = i[:getID]()
            if l in edgeids
                push!(inroads, l2id[l])
            end
        end
        outgos = node[:getOutgoing]()
        outroads = []
        for o in outgos
            l = o[:getID]()
            if l in edgeids
                push!(outroads, l2id[l])
            end
        end
        jtype = node[:getType]()
        clist = node[:getConnections]() # list of all connections
        reds = zeros(length(clist))
        cycle = 10
        if contains(jtype, "traffic_light")
            cycle = 0
            tl = net[:getTLSSecure](nodeid)
            ps = ((tl[:getPrograms]())["0"])[:getPhases]()
            for p in ps
                light = p[1]
                duration = p[2]
                cycle += duration
                for (lidx, l) in enumerate(light)
                    if l == 'r'
                        reds[lidx] += duration
                    elseif l == 's'
                        reds[lidx] += duration*(2/3)
                    end
                end
            end
        elseif contains(jtype, "stop")
            reds = 5*ones(length(clist))
        end

        cdict = Dict()
        for (cidx, c) in enumerate(clist)
            fromEdge = c[:getFrom]()
            toEdge = c[:getTo]()
            if fromEdge[:getID]() in edgeids && toEdge[:getID]() in edgeids
                spdlim = fromEdge[:getSpeed]()
                fromId = l2id[fromEdge[:getID]()]
                toId = l2id[toEdge[:getID]()]
                direction = c[:getDirection]()
                # Saturation equation S = 990 + 288TL + 8.5SL - 26.8G
                if (fromId, toId) in keys(cdict)
                    throughlane = cdict[(fromId, toId)].throughlane + 1
                    pp = (990 + 288*throughlane + 8.5*spdlim*3.6)*throughlane
                    if direction == 'r'
                        pp = pp*2/3
                    elseif direction == 'l'
                        pp = pp*2/3
                    end
                    cdict[(fromId, toId)] = Connection(fromId, toId, throughlane, pp, direction, reds[cidx]) #vph
                else
                    pp = (990 + 288*1 + 8.5*spdlim*3.6) #vph
                    if direction == 'r'
                        pp = pp*2/3
                    elseif direction == 'l'
                        pp = pp*2/3
                    end
                    cdict[(fromId, toId)] = Connection(fromId, toId, 1, pp, direction, reds[cidx])
                end
            end
        end
        j = Junction(index, nodeid, jtype, inroads, outroads, cdict, cycle)
        push!(js, j)
    end
    es = []
    for (index, edgeid) in enumerate(edgeids)
        edge = net[:getEdge](edgeid)
        etype = etypes[index]
        if etype == "major"
            L = edge[:getLength]()
            numlane = length(edge[:getLanes]())
            storage = Int(floor(L * numlane / CARLENGTH))
            tonode = edge[:getToNode]()
            co = Dict{Int, Array{Flow,1}}()
            for outroad in edge[:getOutgoing]()
                outroadid = outroad[1][:getID]()
                if outroadid in edgeids
                    co[l2id[outroadid]] = Flow[]
                end
            end
            spdlim = edge[:getSpeed]()
            e = Urban(index, storage, co, storage, 1, spdlim, L, numlane)
            push!(es, e)
        end
    end
    return Network(es, js)
end
function Base.show(io::IO, m::Network)
    print(io, "A road network")
end
end
