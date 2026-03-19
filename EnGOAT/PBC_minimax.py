import heapq

def minimax_periodic(adj, source, target_offset):

    source_state = (source, 0)

    cost = {source_state: 0.0}

    # parent stores: previous_state, transition_id, boundary_crossing
    parent = {source_state: (None, None, 0)}

    pq = [(0.0, source_state)]

    OFFSET_LIMIT = abs(target_offset) + 2  # safety margin

    while pq:
        c, (u, offset) = heapq.heappop(pq)

        if cost[(u, offset)] < c:
            continue

        # Goal reached
        if u == source and offset == target_offset:

            path = []
            state = (u, offset)

            while state is not None:
                node, off = state
                prev_state, trans_id, crossing = parent[state]

                path.append((node, trans_id, crossing))

                state = prev_state

            path.reverse()
            return c, path

        for v, w, shift, trans_id in adj[u]:

            new_offset = offset + shift

            if abs(new_offset) > OFFSET_LIMIT:
                continue

            state_v = (v, new_offset)

            new_cost = max(c, w)

            if state_v not in cost or new_cost < cost[state_v]:

                cost[state_v] = new_cost

                # shift already represents boundary crossing info
                parent[state_v] = ((u, offset), trans_id, shift)

                heapq.heappush(pq, (new_cost, state_v))

    return float("inf"), None