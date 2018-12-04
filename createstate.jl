function create_state!(net, condemands, queuedict, urate, rng)
    # put routes inside network
    for route in net.routes
        roadIds = route.roads
        for (idx, rid) in enumerate(roadIds)
            if idx == length(roadIds)
                nexto = -1
            else
                nexto = roadIds[idx+1]
            end
            road = get_road(net, rid)
            if nexto in keys(road.collection)
                hasroute = false
                for flow2 in road.collection[nexto]
                    if flow2.routeid == route.id
                        hasroute = true
                        break
                    end
                end
                if !hasroute
                    push!(road.collection[nexto], Flow(routeid = route.id, roadindex = nexto == -1 ? -1:idx+1))
                end
            else
                road.collection[nexto] = [Flow(routeid = route.id, roadindex = nexto == -1 ? -1:idx+1)]
            end
        end
    end
    # first unpack the queues
    for (firstedgeid, queue) in queuedict
        #firstedge is the queuue starter
        firstedge = get_road(net, firstedgeid)
        if queue > firstedge.maxstorage
            realqueue =firstedge.maxstorage
            firstedge.storage = 0.0
            N = 0
            rr = []
            for (nexto,fs) in firstedge.collection
                for f in fs
                    N += 1
                    routeid = f.routeid
                    roadindex = f.roadindex
                    if roadindex != -1
                        push!(rr, (routeid, roadindex-2))
                    else
                        if length(net.routes[routeid].roads) >= 2
                            push!(rr, (routeid, length(net.routes[routeid].roads)-1))
                        end
                    end
                end
            end
            assert(N > 0)
            share = realqueue/N
            for (nexto,fs) in firstedge.collection
                for f in fs
                    f.queue = share
                end
            end
            firstedge.delay = 1
            queue = queue-realqueue

            routesharequeue = queue/length(rr)
            for (rotid, prid) in rr
                previndex = prid
                subqueue = routesharequeue
                curroadid = firstedgeid
                while subqueue != 0.0
                    prevroadid = previndex==0?-1:net.routes[rotid].roads[previndex]
                    if prevroadid == -1
                        update_vehicles!(net, condemands, rotid, net.routes[rotid].roads[1], subqueue)
                        subqueue = 0.0
                        break
                    else
                        proad = get_road(net, prevroadid)
                        if subqueue > proad.maxstorage
                            subqueue = subqueue - proad.maxstorage
                            proad.storage = 0.0
                            realsubqueue = proad.maxstorage
                            for f in proad.collection[curroadid]
                                if f.routeid == rotid
                                    f.queue = realsubqueue
                                    break
                                end
                            end
                            previndex -= 1
                            curroadid = prevroadid
                        else
                            realsubqueue = subqueue
                            proad.storage = proad.maxstorage - realsubqueue
                            for f in proad.collection[curroadid]
                                if f.routeid == rotid
                                    f.queue = realsubqueue
                                    break
                                end
                            end
                            subqueue = 0.0
                            break
                        end

                    end
                end
            end

        else
            realqueue = queue
            firstedge.storage = firstedge.maxstorage - realqueue
            N = 0
            for (nexto,fs) in firstedge.collection
                for f in fs
                    N += 1
                end
            end
            assert(N > 0)
            share = realqueue/N
            for (nexto,fs) in firstedge.collection
                for f in fs
                    f.queue = share
                end
            end
        end
    end
    # update arrival2
    for route in net.routes
        q = Inf
        roadIds = route.roads
        for (idx, rid) in enumerate(roadIds)
            road = get_road(net, rid)
            if idx == length(roadIds)
                nexto = -1
            else
                nexto = roadIds[idx+1]
            end

            flow = nothing
            for (flowidx, f) in enumerate(road.collection[nexto])
                if f.routeid == route.id
                    flow = f
                    generate_arrival2!(road, nexto, flowidx, q)
                    break
                end
            end
            if flow.queue > 0.0
                j = net.junctions[road.D]
                qq = get_capacity(j, rid, nexto)
                if qq < q
                    q = qq
                end
            end
        end
    end
    # update uncon-user demand
    for rate in keys(urate)
        route1,route2,fastest=urate[rate]
        time1 = compute_time(net, condemands, route1)
        time2 = compute_time(net, condemands, route2)
        dist1 = compute_distance(net, route1)
        dist2 = compute_distance(net, route2)
        percent1 = 1.0
        percent2 = 0.0
        if dist1 <= dist2
            if time1 <= time2
                percent1 = 1.0
                percent2 = 0.0
            else
                percent1 = 1-fastest
                percent2 = fastest
            end
        else
            if time1 <= time2
                percent1 = fastest
                percent2 = 1-fastest
            else
                percent1 = 0.0
                percent2 = 1.0
            end
        end
        condemands[(net.routes[route1].O, net.routes[route1].roads[1])][1].rate += rate*percent1
        condemands[(net.routes[route2].O, net.routes[route2].roads[1])][1].rate += rate*percent2
    end
    state1 = deepcopy(net)
    return state1
end
function state_step!(net, condemands, queuedict, urate, rng)
    rate1 = condemands[(net.routes[1].O,net.routes[1].roads[1])][1].rate
    rate2 = condemands[(net.routes[2].O,net.routes[2].roads[1])][1].rate
    state1 = create_state!(net, condemands, queuedict, urate, rng)
    #time1_1 = compute_time(state1, condemands, state1.routes[1].id)
    #time1_2 = compute_time(state1, condemands, state1.routes[2].id)

    energy1_1 = compute_energy_cost(state1, state1.routes[1].id)
    energy1_2 = compute_energy_cost(state1, state1.routes[2].id)

    for i = 1:150
        step(net, condemands, rng)
    end
    state2 = deepcopy(net)
    energy2_1 = compute_energy_cost(state2, state2.routes[1].id)
    energy2_2 = compute_energy_cost(state2, state2.routes[2].id)


    return state2, (energy1_1 + energy2_1)*0.5*rate1 + (energy1_2 + energy2_2)*0.5*rate2
end
function update_vehicles!(net, demands, rotid, firstedgeid, subqueue)
    route = get_route(net, rotid)
    O = route.O
    for d in demands[(O, firstedgeid)]
        if d.routeid == rotid
            d.vehicles += subqueue
            break
        end
    end
end
function generate_arrival2!(road, nexto, flowidx, q)
    if q == Inf
        qq = 0.0
        if typeof(road) == Freeway
            road.delay = Int(ceil(road.storage*1.5*CARLENGTH/(road.spdlim*Tu*(road.numlane/4+0.75))))
        else
            road.delay = Int(ceil(road.storage*CARLENGTH/(road.spdlim*Tu*road.numlane)))
        end
    else
        if typeof(road) == Freeway
            vdelay = Int(ceil(road.storage*1.5*CARLENGTH/(road.spdlim*Tu*(road.numlane/4+0.75))))
        else
            vdelay = Int(ceil(road.storage*CARLENGTH/(road.spdlim*Tu*road.numlane)))
        end
        while road.storage > 0 && vdelay >= 1
            road.collection[nexto][flowidx].arrival2[vdelay] = road.storage > q*Tu ? q*Tu : road.storage
            road.storage = max(0.0, road.storage-q*Tu)
            vdelay -= 1
        end
        if typeof(road) == Freeway
            road.delay = Int(ceil(road.storage*1.5*CARLENGTH/(road.spdlim*Tu*(road.numlane/4+0.75))))
        else
            road.delay = Int(ceil(road.storage*CARLENGTH/(road.spdlim*Tu*road.numlane)))
        end
    end
end
