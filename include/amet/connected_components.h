#pragma once

#include <graph/graph.hpp>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <limits>

namespace amet {

    /**
     * @brief Find the weakly/strongly connected components of an adjacency graph.
     * 
     * @tparam G               The concept-constrained adjacency list graph.
     * @tparam ComponentMap    A standard indexable lookup container (std::map or std::vector).
     * 
     * @return std::size_t     The total number of isolated component networks found.
     */
    template <graph::adjacency_list G, typename ComponentMap>
    std::size_t connected_components(G&& g, ComponentMap& component_map) {
        using vid_t = graph::vertex_id_t<G>;
        
        const std::size_t N = graph::num_vertices(g);
        if (N == 0) return 0;

        constexpr std::size_t UNVISITED = std::numeric_limits<std::size_t>::max();
        
        // 1. Initialize all active nodes as unvisited inside your lookup
        for (auto&& [uid, u] : graph::views::vertexlist(g)) {
            component_map[uid] = UNVISITED;
        }

        std::size_t component_id = 0;
        std::queue<vid_t> q;

        // 2. Loop over every vertex to capture disconnected or isolated clusters cleanly
        for (auto&& [start_id, start_vertex] : graph::views::vertexlist(g)) {
            if (component_map[start_id] != UNVISITED) continue;

            // Found a new isolated root cluster!
            component_map[start_id] = component_id;
            q.push(start_id);

            while (!q.empty()) {
                vid_t curr_id = q.front();
                q.pop();

                // Locate the active vertex token safely via the core CPO lookup
                auto&& curr_vertex = *graph::find_vertex(g, curr_id);

                // Traverse neighbors explicitly via top-level CPOs
                for (auto&& edge_ref : graph::views::incidence(g, curr_vertex)) {
                    vid_t neighbor_id = graph::target_id(g, edge_ref);

                    if (component_map[neighbor_id] == UNVISITED) {
                        component_map[neighbor_id] = component_id;
                        q.push(neighbor_id);
                    }
                }
            }
            
            // Advance the allocation key for the next isolated structural block
            ++component_id;
        }

        return component_id; // Returns total component count natively
    }
}
