include("constants.jl")

function compute_default_time(net, demands, routeid)
    route = get_route(net, routeid)
    roadids = route.roads
    traveltime = 0.0
    for (idx, roadId) in enumerate(roadids)
        road = get_road(net, roadId)
        traveltime += road.length/road.spdlim
    end
    return traveltime
end
function compute_time(net, demands, routeid)
    route = get_route(net, routeid)
    origin = route.O
    roadids = route.roads
    outid =  roadids[1]
    ds = get_source(demands, origin, outid)
    demand = nothing
    for d in ds
        if d.routeid == routeid
            demand = d
            break
        end
    end
    traveltime = 0.0
    curqueue = 0.0
    curqueue += demand.vehicles
    curQ = demand.rate
    for (idx, roadId) in enumerate(roadids)
        road = get_road(net, roadId)
        nexto = -1
        if idx < length(roadids)
            nexto = roadids[idx+1]
        end
        flows = []
        if nexto in keys(road.collection)
            flows = road.collection[nexto]
        end
        delay = road.delay*Tu
        thisqueue = 0
        for f in flows
            thisqueue += f.queue
        end
        thisQ = get_capacity(net.junctions[road.D], roadId, nexto)
        if road.storage <= 1 # this means over crowded
            #v = 1/((1/road.spdlim)+(ra/(CARLENGTH*satflow*(cycle-ra))))
            curQ = thisQ
            curqueue += thisqueue

        else
            traveltime +=  curqueue/curQ
            traveltime += delay
            curqueue = thisqueue
            curQ = thisQ
        end
    end
    traveltime += curqueue/curQ
    return traveltime
end
function compute_distance(net, routeid)
    route = get_route(net, routeid)
    roadids = route.roads
    traveldistance = 0.0
    for (idx, roadId) in enumerate(roadids)
        road = get_road(net, roadId)
        traveldistance += road.length
    end
    return traveldistance
end
function compute_energy_cost(net, routeid)
    route = get_route(net, routeid)
    roadids = route.roads
    travelenergy = 0.0
    curtraveltime = 0.0
    curqueue = 0.0
    curqueuedistance = 0.0
    curQ = 0.5
    curV = 17.88
    curspdlim = 17.88
    for (idx, roadId) in enumerate(roadids)
        road = get_road(net, roadId)
        numlane= road.numlane
        spdlim = road.spdlim
        nexto = -1
        if idx < length(roadids)
            nexto = roadids[idx+1]
        end
        flows = []
        if nexto in keys(road.collection)
            flows = road.collection[nexto]
        end
        delay = road.delay*Tu
        thisqueuedistance = 0.0
        thisqueue = 0.0
        for f in flows
            thisqueue += f.queue
            if typeof(road) == Freeway
                thisqueuedistance += f.queue*CARLENGTH*1.5/(numlane/4 + 0.75)
            else
                thisqueuedistance += f.queue*CARLENGTH/numlane
            end
        end
        thisQ = get_capacity(net.junctions[road.D], roadId, nexto)
        #thisV = get_velocity(net.junctions[road.D], roadId, nexto)
        if road.storage <= 1.0 # this means over crowded
            #v = 1/((1/road.spdlim)+(ra/(CARLENGTH*satflow*(cycle-ra))))
            curQ = thisQ
            curqueuedistance += thisqueuedistance
            curqueue += thisqueue
            curspdlim = spdlim
            if curqueue == 0.0
                curV = (curqueuedistance)/((curqueue/curQ) + 1e-5)
            else
                curV = (curqueuedistance)/(curqueue/curQ)
            end
            thisV = get_velocity(net, net.junctions[road.D], roadId, nexto)
            #println("Vel of ", roadId, " to ", nexto, " is: ", curV, ".Or Wardrop: ", thisV, " Spdlim is: ", curspdlim)
            #println("QDistance: ", curqueuedistance, " FDistance: ",road.length-thisqueuedistance,  " Total distance: ", road.length, " FErate: " , compute_energy_rate(curspdlim, curspdlim), " Erate: ", compute_energy_rate(curV, curspdlim))
        else
            travelenergy +=  curqueuedistance*compute_energy_rate(curV, curspdlim)
            travelenergy += max(road.length-thisqueuedistance, 0)*compute_energy_rate(curspdlim, curspdlim)
            curQ = thisQ
            curqueuedistance = thisqueuedistance
            curqueue = thisqueue
            curspdlim = spdlim
            if curqueue == 0.0
                curV = (curqueuedistance)/((curqueue/curQ) + 1e-5)
            else
                curV = (curqueuedistance)/(curqueue/curQ)
            end
            thisV = get_velocity(net, net.junctions[road.D], roadId, nexto)
            #println("Vel of ", roadId, " to ", nexto, " is: ", curV, ".Or Wardrop: ", thisV, " Spdlim is: ", curspdlim)
            #println("QDistance: ", curqueuedistance, " FDistance: ",road.length-thisqueuedistance,  " Total distance: ", road.length, " FErate: " , compute_energy_rate(curspdlim, curspdlim), " Erate: ", compute_energy_rate(curV, curspdlim))

        end
    end
    curV = 0.0
    if curqueue == 0.0
        curV = (curqueuedistance)/((curqueue/curQ) + 1e-5)
    else
        curV = (curqueuedistance)/(curqueue/curQ)
    end
    travelenergy += curqueuedistance*compute_energy_rate(curV, curspdlim)
    return travelenergy
end
function get_velocity(net, junction, oid, did)
    oroad = get_road(net, oid)
    if did == -1
        return oroad.spdlim
    end
    spdlim = oroad.spdlim
    droad = get_road(net, did)
    numlane = oroad.numlane
    k1 = numlane/CARLENGTH
    if typeof(oroad) == Ramp && typeof(droad) == Freeway
        spdlim = 29.06
    elseif typeof(oroad) == Freeway && typeof(droad) == Ramp
        spdlim = 17.88
        k1 = (numlane/4+0.75)/(CARLENGTH*1.5)
    end
    con = junction.connections[(oid, did)]
    p = con.satflow # veh/hr
    p = p/3600
    r = con.red
    c = junction.cycle
    ra = r+10 >= c ? r : r+10
    #println(p, " ", ra, " ", k1, " ", (k1*ra)/(p*(c-ra)) )
    v0 = 1/(1/spdlim + (k1*ra)/(p*(c-ra)))
    return v0
end
function compute_energy_rate(v, spdlim)
    return EVitp(spdlim, v)
end
