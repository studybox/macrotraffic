module discretemacrotraffic

using PyCall
using DataStructures
using Distributions
using StaticArrays
using Parameters
using Graphs
import Base: show, step, length
export Flow,
       Urban,
       Freeway,
       Ramp,
       Junction,
       Connection,
       Demand,
       Route,
       Network,
       import_network,
       add_route!,
       step,
       generate_demands
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
    O::Int
    D::Int
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

mutable struct Freeway
    id::Int
    O::Int
    D::Int
    storage::Int
    collection::Dict{Int, Array{Flow,1}}
    maxstorage::Int
    delay::Int
    spdlim::Float64
    length::Float64
    numlane::Int
end
mutable struct Ramp
    id::Int
    O::Int
    D::Int
    storage::Int
    collection::Dict{Int, Array{Flow,1}}
    maxstorage::Int
    delay::Int
    spdlim::Float64
    length::Float64
    numlane::Int
end

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
    rate::Float64
    O::Int
    D::Int
    routeid::Int
    vehicles::Int = 0
end




const Road = Union{Urban, Freeway, Ramp}

function get_source(demands::DefaultDict{Tuple{Int, Int}, Array{Demand,1}}, junctionid::Int, outroadid::Int)
    return demands[(junctionid, outroadid)]
end

struct Network
    roads::SVector
    junctions::SVector
    n2id::Dict{String, Int}
    l2id::Dict{String, Int}
    routes::Array{Route, 1}
end
Network(rs::Array, js::Array, n2id::Dict, l2id::Dict) = Network(SVector{length(rs)}(rs), SVector{length(js)}(js), n2id, l2id, Array{Int, 1}[])
function get_road(net::Network, id::Int)
    return net.roads[id]
end
function get_flow1(road, did, Qcap)
    flows = road.collection[did]
    q = 0
    for f in flows
        q += f.queue + f.arrival2[0]
    end
    return min(q, Int(ceil(Tu*Qcap)))
end
function get_capacity(junction::Junction, oid::Int, did::Int)
    con = junction.connections[(oid, did)]
    p = con.satflow # veh/hr
    r = con.red
    c = junction.cycle
    ra = r+10 >= c ? r : r+10
    Q = p*(c-ra)/c/3600
    #println("Cap of ", oid, " to ", did, " is: ", Tu*Q)

    return Q
end
function create_vehicles!(demands::Array{Demand, 1}, rng::AbstractRNG)
    pq = PriorityQueue{Int, Int}()
    totalf = 0
    for (idx, d) in enumerate(demands)
        f = create_vehicles!(d, rng)
        enqueue!(pq, idx=>f)
        totalf += f
    end
    return (pq, totalf)
end
function create_vehicles!(demand::Demand, rng::AbstractRNG)
    u = Poisson(Tu*demand.rate)
    newveh = rand(u)
    demand.vehicles += newveh
    return demand.vehicles
end

function update_storage!(road::Road, fs::PriorityQueue, tf::Int)

    if tf <= road.storage
        road.storage -= tf
        return fs
    else
        pq = PriorityQueue{Int, Int}()
        N = length(fs)
        for i = 1:N
            (idx, f) = peek(fs)
            if f <= road.storage/(N-i+1)
                road.storage -= f
                enqueue!(pq, idx=>f)
            else
                #println("storage: ", road.storage)
                #println("N: ", N-i+1)
                res = road.storage % (N-i+1)
                #println("res: ", res)
                share = Int(floor(road.storage / (N-i+1)))
                #println("share: ", share)
                fff = 0
                for (ii, idx2) in enumerate(keys(fs))
                    q = ii > (N-i+1)-res ? share+1 : share
                    fff+=q
                    enqueue!(pq, idx2=>q)
                end
                #println("No more room for demands")

                assert(fff == road.storage)
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
        if N == 0
            break
        end
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
length(route::Route) = length(route.roads)
function update_arrival1!(road, route, routeid, roadindex, arrival1)
    if length(route) == roadindex
        nexto = - 1
    else
        nexto = route.roads[roadindex+1]
    end
    #println(nexto)
    if nexto in keys(road.collection)
        hasroute = false

        for flow2 in road.collection[nexto]
            if flow2.routeid == routeid
                flow2.arrival1 += arrival1
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

function get_route(net, routeid)
    return net.routes[routeid]
end
function update_arrival1!(road, demand, ff, net)
    for (idx, d) in enumerate(demand)
        arrival1 = ff[idx]
        routeid = d.routeid
        roadindex = 1
        route = get_route(net, routeid)
        update_arrival1!(road, route, routeid, roadindex, arrival1)
    end
    #println("arrival1 : ", road.collection)
end

function update_storage!(road)
    for d in keys(road.collection)
        if d == -1
            for f in road.collection[d]
                road.storage += f.arrival2[0]
            end
        else
            for f in road.collection[d]
                road.storage += f.depart
            end
        end
    end
    road.storage = min(road.storage, road.maxstorage)
end

function update_queue!(road)
    for d in keys(road.collection)
        if d != -1
            for f in road.collection[d]
                f.queue = max(f.queue + f.arrival2[0] - f.depart , 0)
            end
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
                if t-1 >= 0
                    new_arrival2[t-1] = f.arrival2[t]
                end
            end
            f.arrival2 = new_arrival2
        end
    end
    if typeof(road) == Freeway
        road.delay = Int(ceil(road.storage*1.5*CARLENGTH/(road.spdlim*Tu*(road.numlane/4+0.75))))
    else
        road.delay = Int(ceil(road.storage*CARLENGTH/(road.spdlim*Tu*road.numlane)))
    end
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
                #println("flows2 in demands: ", flows1)
                # Must satisfy demand first
                # Compute and update remaining storage
                flows2 = update_storage!(outroad, flows1, totalflows1)
                #println("flows2 in act demands: ", flows2)
                update_arrival1!(outroad, demand, flows2, network)
                # update demand if not all satisfied
                update_vehicles!(demand, flows2)
            end
            flows1 = PriorityQueue{Int, Int}()
            totalflows1 = 0
            for inid in junction.inroads
                if is_connected(junction, inid, outid)
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
function add_route!(net::Network, route::Array{String,1})
    O = net.roads[net.l2id[route[1]]].O
    D = net.roads[net.l2id[route[end]]].D
    ll = []
    for r in route
        l = net.l2id[r]
        push!(ll, l)
    end
    id = length(net.routes)+1
    route = Route(id, ll, O, D)
    add_route!(net, route)
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
                if is_connected(junction, edgeidx, outid)
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
                    #=
                    if index == 23
                        println("Nodeid: ", nodeid)
                        println("lights: ", light)
                        println("l: ", l)
                        println("is red?: ", l == 'r')
                    end
                    =#
                    if l == 'r'
                        reds[lidx] += duration
                    elseif l == 's'
                        reds[lidx] += duration*(2/3)
                    end
                end
            end

        elseif contains(jtype, "allway_stop")
            reds = 5*ones(length(clist))
        end

        cdict = Dict()
        fltdset = Set()
        for (cidx, c) in enumerate(clist)
            fromEdge = c[:getFrom]()
            fromLane = c[:getFromLane]()
            toEdge = c[:getTo]()
            stateid = c[:getTLLinkIndex]() + 1
            if stateid <= 0
                stateid = cidx
            end
            if fromEdge[:getID]() in edgeids && toEdge[:getID]() in edgeids
                spdlim = fromEdge[:getSpeed]()

                fromId = l2id[fromEdge[:getID]()]
                fromLaneId = fromLane[:getID]()
                toId = l2id[toEdge[:getID]()]
                if etypes[fromId] == "ramp" && etypes[toId] == "expressway"
                    spdlim = 29.06
                end

                if etypes[fromId] == "expressway" && etypes[toId] == "ramp"
                    spdlim = 17.88
                end

                direction = c[:getDirection]()
                fltd = (fromLaneId,direction)
                # Saturation equation S = 990 + 288TL + 8.5SL - 26.8G
                if (fromId, toId) in keys(cdict)
                    if !(fltd in fltdset)
                        throughlane = cdict[(fromId, toId)].throughlane + 1
                        pp = (990 + 288*throughlane + 8.5*spdlim*3.6)*throughlane
                        if direction == 'r'
                            pp = pp*2/3
                        elseif direction == 'l'
                            pp = pp*2/3
                        end

                        cdict[(fromId, toId)] = Connection(fromId, toId, throughlane, pp, direction, reds[stateid]) #vph
                    end
                else
                    push!(fltdset, fltd)
                    pp = (990 + 288*1 + 8.5*spdlim*3.6) #vph
                    if direction == 'r'
                        pp = pp*2/3
                    elseif direction == 'l'
                        pp = pp*2/3
                    end
                    if (fromId, toId) == (46, 47)
                        println(cidx)
                        println(reds[cidx])
                    end
                    cdict[(fromId, toId)] = Connection(fromId, toId, 1, pp, direction, reds[stateid])
                end
            end
        end
        j = Junction(index, nodeid, jtype, inroads, outroads, cdict, cycle)
        push!(js, j)
    end
    es = []
    for (index, edgeid) in enumerate(edgeids)
        etype = etypes[index]
        if etype == "MAJOR"
            edge = net[:getEdge](edgeid)
            L = edge[:getLength]()
            Oname = edge[:getFromNode]()[:getID]()
            Dname = edge[:getToNode]()[:getID]()
            O = n2id[Oname]
            D = n2id[Dname]
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
            delay = Int(ceil(storage*CARLENGTH/(spdlim*Tu*numlane)))
            e = Urban(index, O, D, storage, co, storage, delay, spdlim, L, numlane)
            push!(es, e)
        elseif etype == "EXPRESSWAY"
            edge = net[:getEdge](edgeid)
            L = edge[:getLength]()
            Oname = edge[:getFromNode]()[:getID]()
            Dname = edge[:getToNode]()[:getID]()
            O = n2id[Oname]
            D = n2id[Dname]
            numlane = length(edge[:getLanes]())

            storage = Int(floor(L/ (CARLENGTH*1.5))) + Int(floor(L*(numlane-1)*0.25/(CARLENGTH*1.5)))
            tonode = edge[:getToNode]()
            co = Dict{Int, Array{Flow,1}}()
            for outroad in edge[:getOutgoing]()
                outroadid = outroad[1][:getID]()
                if outroadid in edgeids
                    co[l2id[outroadid]] = Flow[]
                end
            end
            spdlim = edge[:getSpeed]()
            delay = Int(ceil(L/(spdlim*Tu)))
            e = Freeway(index, O, D, storage, co, storage, delay, spdlim, L, numlane)
            push!(es, e)
        elseif etype == "RAMP" || etype == "COLLECTOR"
            edge = net[:getEdge](edgeid)
            L = edge[:getLength]()
            Oname = edge[:getFromNode]()[:getID]()
            Dname = edge[:getToNode]()[:getID]()
            O = n2id[Oname]
            D = n2id[Dname]
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
            delay = Int(ceil(storage*CARLENGTH/(spdlim*Tu*numlane)))
            e = Ramp(index, O, D, storage, co, storage, delay, spdlim, L, numlane)
            push!(es, e)
        end
    end
    return Network(es, js, n2id, l2id)
end
function Base.show(io::IO, r::Road)
    println(io, "id: ", r.id)
    println(io, "O-D: ", r.O, " ", r.D)
    println(io, "storage: ", r.storage, "veh")
    println(io, "Time to Queue: ", r.delay*Tu, "s")
    for (outrid, flows) in r.collection
        println(io, "Flows to road ", outrid, ":")
        for f in flows
            Base.show(io, f)
        end
    end
end
function Base.show(io::IO, f::Flow)
    println(io, "    routeid: ", f.routeid)
    println(io, "    queue: ", f.queue, "veh")
    println(io, "    arrival: ", f.arrival2[0], "veh")
    aa = 0
    for (t, a) in f.arrival2
        aa += a
    end
    println(io, "    future arrivals: ", aa, "veh")
    println(io, "    departure: ", f.depart, "veh")
end
function Base.show(io::IO, m::Network)
    print(io, "A road network")
end
end
#=
edgesid = edgeid
LL = []
numlane = []
storage = 0
co = Dict{Int, Array{Flow,1}}()
for (eidx, eid) in enumerate(edgesid)
    edge = net[:getEdge](eid)
    L  = edge[:getLength]()
    nl = length(edge[:getLanes]())
    push!(LL, L)
    push!(numlane, nl)
    if eidx == 1
        Oname = edge[:getFromNode]()[:getID]()
        O = n2id[Oname]
        storage += Int(floor(L/CARLENGTH))+Int(floor(L*(numlane-1)*0.5/CARLENGTH))
        spdlim = edge[:getSpeed]()
    elseif eidx == length(edgesid)
        Dname = edge[:getToNode]()[:getID]()
        D = n2id[Dname]
        storage += Int(floor(L/ CARLENGTH))
        tonode = edge[:getToNode]()
        for outroad in edge[:getOutgoing]()
            outroadid = outroad[1][:getID]()
            if outroadid in edgeids
                co[l2id[outroadid]] = Flow[]
            end
        end
    else
        storage += Int(floor(L/CARLENGTH))+Int(floor(L*(numlane-1)*0.5/CARLENGTH))
    end
end
delay = Int(ceil(storage*CARLENGTH/(spdlim*Tu*numlane)))
e = Freeway(index, O, D, storage, co, storage, delay, spdlim, LL, numlane)
push!(es, e)
=#
