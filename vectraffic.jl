demand, depart, arrival1, arrival2, storage, queue = initilize(S)
for idx = 1:H
    update!(demand, A)
    arrival1 = arrival1 + DE_A1*demand
    depart = get_outgoing_flow(arrival2, queue)
    real_depart!(depart, storage)
    arrival1 = real_arrival1(depart)
    update_storage!(storage)
    update_queue!(queue)
    update_arrival2!(arrival2)
end


update!(demand, A)
