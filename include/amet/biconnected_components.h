#pragma once

#include <graph/graph.hpp> // Bring in only the standard concept definitions
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <algorithm>
#include <limits>

namespace amet {

    // A clean, explicitly typed CPO-styled implementation
    template <graph::adjacency_list G, typename OuterContainer>
    void biconnected_components(G&& g, OuterContainer& components) {
        using vid_t = graph::vertex_id_t<G>;
        using inner_type = typename OuterContainer::value_type;
        
        const std::size_t N = graph::num_vertices(g);
        if (N == 0) return;

        constexpr std::size_t UNVISITED = std::numeric_limits<std::size_t>::max();

        // Use clean, explicit local tracking vectors instead of library macro property maps
        std::map<vid_t, std::size_t> disc;
        std::map<vid_t, std::size_t> low;
        std::map<vid_t, vid_t> parent;
        
        for (auto&& [uid, u] : graph::views::vertexlist(g)) {
            disc[uid] = UNVISITED;
            low[uid] = UNVISITED;
        }

        std::size_t timer = 0;
        std::stack<std::pair<vid_t, vid_t>> edge_stack;

        // DFS Frame tracking state explicitly without hiding under anonymous auto variables
        struct dfs_frame {
            vid_t uid;
            decltype(graph::views::incidence(g, std::declval<graph::vertex_t<G>>())) edge_view;
            typename decltype(graph::views::incidence(g, std::declval<graph::vertex_t<G>>()))::iterator it;
        };
        
        std::stack<dfs_frame> dfs_stack;

        // Helper lambda to flush components cleanly using explicit standard types
        auto flush_component = [&](vid_t u, vid_t v) {
            std::set<vid_t> unique_vertices;
            while (!edge_stack.empty()) {
                auto [eu, ev] = edge_stack.top();
                edge_stack.pop();
                unique_vertices.insert(eu);
                unique_vertices.insert(ev);
                if (eu == u && ev == v) break;
            }
            
            inner_type comp;
            comp.reserve(unique_vertices.size());
            for (const vid_t vv : unique_vertices) {
                comp.push_back(static_cast<typename inner_type::value_type>(vv));
            }
            components.push_back(std::move(comp));
        };

        // Outer loop to naturally process disconnected graphs or periodic boundary sets
        for (auto&& [start_id, start_vertex] : graph::views::vertexlist(g)) {
            if (disc[start_id] != UNVISITED) continue;

            disc[start_id] = low[start_id] = timer++;
            auto start_edges = graph::views::incidence(g, start_vertex);
            
            dfs_stack.push(dfs_frame{
                .uid = start_id,
                .edge_view = start_edges,
                .it = start_edges.begin()
            });

            while (!dfs_stack.empty()) {
                auto& frame = dfs_stack.top();
                auto uid = frame.uid;
                
                if (frame.it == frame.edge_view.end()) {
                    // Backtrack step
                    dfs_stack.pop();
                    if (!dfs_stack.empty()) {
                        auto& parent_frame = dfs_stack.top();
                        vid_t pid = parent_frame.uid;
                        
                        low[pid] = std::min(low[pid], low[uid]);
                        
                        if (low[uid] >= disc[pid]) {
                            flush_component(pid, uid);
                        }
                    }
                    continue;
                }

                // UNIFYING STEP: Access the neighbor data explicitly via standard CPOs
                auto&& edge_ref = *(frame.it);
                vid_t vid = graph::target_id(g, edge_ref);
                ++(frame.it); // Advance frame iterator explicitly

                if (vid == uid) continue; // Skip self-loops safely

                if (disc[vid] == UNVISITED) {
                    parent[vid] = uid;
                    disc[vid] = low[vid] = timer++;
                    edge_stack.push({uid, vid});
                    
                    auto next_vertex = *graph::find_vertex(g, vid);
                    auto vid_edges = graph::views::incidence(g, next_vertex);
                    
                    dfs_stack.push(dfs_frame{
                        .uid = vid,
                        .edge_view = vid_edges,
                        .it = vid_edges.begin()
                    });
                } else if (parent.contains(uid) && parent[uid] == vid) {
                    continue; // Skip direct link back to DFS parent
                } else {
                    // True back-edge logic
                    if (disc[vid] < disc[uid]) {
                        edge_stack.push({uid, vid});
                    }
                    low[uid] = std::min(low[uid], disc[vid]);
                }
            }
        }
    }
}
