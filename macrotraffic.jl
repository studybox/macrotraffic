module macrotraffic

using Distributions
using StaticArrays
using Parameters
using Graphs
#import
#export

struct Flow{T}
    id::Int
    vehicles::Int
    route::Array{T, 1}
end

struct Urban
    id::Int
    queue::Array{Int, 1}
    storage::Int
end
struct Freeway
end
struct Ramp
end
struct Junction
    id::Int
    connections::Array{,1}
end

@with_kw struct Demand
    rate::Int = -1
    O::Int = -1
    D::Int = -1
    routeid::Int = -1
    vehicles::Int = -1
end
void_demands = Array{Demand, 1}()

function get_source(demands::Dict{Int, Array{Demand,1}}, junction::Junction)
    return demands[Junction.id]
end

const Road = Union{Urban, Freeway, Ramp}

struct Network
    roads::Array{Road, 1}
    junctions::Array{Junction, 1}
end
function get_road(network, id)
end
function get_flow1(road, did, Qcap)
    return minimum(road.queue[did] + road.arrive[did], Tu*Qcap)
end
function get_capacity(junction, oid, did)
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
function step(network, demands, rng)
    for junction in network.junctions
        for outid in junction.outroads
            outroad = get_road(network, outid)
            demand = get_source(demands, junction, outid)
            if length(demand) == 0
                continue
            end
            # compute intentions
            flows1, totalflows1 = create_vehicles!(demand, rng)
            # Must satisfy demand first
            # Compute and update remaining storage
            flows2 = update_storage!(outroad, flows1, totalflows1)
            # update demand if not all satisfied
            update_vehicles!(demand, flows2)
            pq = Priorityqueue{Int, Int}()
            for inid in junction.inroads
                if is_connected(junction.connections, inid, outid)
                    inroad = get_road(network, inid)
                    Qcap = get_capacity(junction, inid, outid)
                    flows1 = get_flow1(inroad, outid, Qcap)
                    enqueue!(pq, inid)
                end
            end
            # TODO decide actual flow
            # TODO update next road
        end
    end
    for road in network.roads
        # TODO Update storage
        # TODO accumulate arrival at link
        # TODO update arrival at queue
        # TODO update queue
    end
end

end
