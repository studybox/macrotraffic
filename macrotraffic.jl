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
    connections::Array{,1}
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
function step(network, demands)
    for junction in network.junctions
        demand = get_source(demands, junction)
        for outid in junction.outroads
            # compute intentions
            flow1 = get_flow1(demand, outid)
            # Must satisfy demand first
            # TODO Compute and update remaining storage
            # TODO update demand if not all satisfied
            for inid in junction.inroads
                if is_connected(junction.connections, inid, outid)
                    inroad = get_road(network, inid)
                    Qcap = get_capacity(junction, inid, outid)
                    flow1 = get_flow1(inroad, outid, Qcap)
                    # TODO save in some priorityqueue for later use
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
