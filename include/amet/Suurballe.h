#pragma once

#include <graph/graph.hpp>
#include <graph/graph_concepts.hpp>
#include <graph/algorithm/dijkstra_shortest_paths.hpp>
#include <vector>
#include <ranges>
#include <optional>
#include <algorithm>

namespace amet{
    template <graph::adj_list::adjacency_list G>
    struct suurballe_result {
        // using edge_t = typename std::remove_cvref_t<G>::vertex_id_type;
        std::vector<graph::adj_list::edge_t<const G>> path1;
        std::vector<graph::adj_list::edge_t<const G>> path2;
        bool success = false;
    };

    template <graph::adj_list::adjacency_list G>
    class Suurballe{
        public:
        Suurballe() = default;
        virtual ~Suurballe() = default;
        template <typename WeightProjection>
        auto operator()(const G& g,
                         graph::vertex_id_t<G> source,
                         graph::vertex_id_t<G> target,
                         WeightProjection weight_proj) -> suurballe_result<G> {
            using vertex_t = typename std::remove_cvref_t<G>::vertex_id_type;
            using edge_t = std::ranges::iterator_t<decltype(graph::adj_list::edges(g, std::declval<const vertex_t&>()))>;
            
            size_t num_v = graph::num_vertices(g);
            suurballe_result<G> result;

            // --- STEP 1: First Dijkstra Pass ---
            std::vector<double> dist1(num_v, std::numeric_limits<double>::max());
            std::vector<graph::vertex_id_t<G>> pred1(num_v);
            
            dist1[source] = 0.0;
            
            // Using standard graph-v3 dijkstra configuration
            // Note: Adjust the parameter shape to match your specific graph-v3 version's CPO signature
            graph::dijkstra_shortest_paths(g, source, graph::container_value_fn(dist1), graph::container_value_fn(pred1), weight_proj);

            // If target is unreachable on the first pass, no disjoint paths exist
            if (dist1[target] == std::numeric_limits<double>::max()) {
                return result;
            }

            // Trace path 1 back to collect its edges
            using edge_nct = edge_t;//std::ranges::iterator_t<decltype(graph::adj_list::edges(const_cast<G&>(g), std::declval<const vertex_t&>()))>;
            std::vector<edge_nct> p1_edges;
            vertex_t curr = target;
            while (curr != source) {
                if (pred1[curr]==pred1[source]) break;
                auto edges_view = graph::adj_list::edges(g, pred1[curr]);
                auto it = std::ranges::find_if(edges_view, [curr, &g](const auto& edge_desc) {
                    return graph::target_id(g, edge_desc) == curr; 
                });
                
                if (it != std::ranges::end(edges_view)) {
                    // Store the actual edge descriptor view/object
                    p1_edges.push_back(it);
                }
                
                curr = pred1[curr];
            }
            std::ranges::reverse(p1_edges);

            // --- STEP 2 & 3: Weight Transformation & Edge Modification View ---
            // Instead of instantiating a completely new modified graph type, we can model
            // the weight updates and edge-reversals through a custom weight function wrapper
            // for the second Dijkstra pass.
            
            auto modified_weight = [&](const edge_t& edge) -> double {
                vertex_t u = graph::source_id(g, *edge);
                vertex_t v = graph::target_id(g, *edge);
                
                // If this edge was used in Path 1, conceptually its reversed counterpart 
                // has a cost of 0, but its forward orientation is blocked / infinite.
                bool is_p1_edge = std::ranges::any_of(p1_edges, [&](const edge_nct& p1_e) {
                    return graph::source_id(g, *p1_e) == u && graph::target_id(g, *p1_e) == v;
                });

                if (is_p1_edge) {
                    return std::numeric_limits<double>::max(); // Block forward path reuse
                }

                double w = static_cast<double>(std::invoke(weight_proj, g, *edge));
                return w - dist1[v] + dist1[u];
            };

            // --- STEP 4: Second Dijkstra Pass ---
            std::vector<double> dist2(num_v, std::numeric_limits<double>::max());
            std::vector<graph::vertex_id_t<G>> pred2(num_v);
            
            dist2[source] = 0.0;
            
            graph::dijkstra_shortest_paths(g, source, graph::container_value_fn(dist2), graph::container_value_fn(pred2), weight_proj);

            if (dist2[target] == std::numeric_limits<double>::max()) {
                // No second edge-disjoint path available
                result.path1.reserve(p1_edges.size());
                for (auto& it : p1_edges) {
                    using target_descriptor_t = typename decltype(result.path1)::value_type;
                    auto const& view = *it;
                    result.path1.push_back(view);
                }
                result.success = false;
                return result;
            }

            // Trace path 2
            std::vector<edge_nct> p2_edges;
            curr = target;
            while (curr != source) {
                if (pred2[curr]==pred2[source]) break;
                    auto edges_view = graph::adj_list::edges(g, pred2[curr]);
                auto it = std::ranges::find_if(edges_view, [curr, &g](const auto& edge_desc) {
                    return graph::target_id(g, edge_desc) == curr; 
                });
                
                if (it != std::ranges::end(edges_view)) {
                    // Store the actual edge descriptor view/object
                    p2_edges.push_back(it);
                }
                
                curr = pred2[curr];
            }
            std::ranges::reverse(p2_edges);

            // TODO: Perform path interlacing (canceling out mutually reversed edges) 
            // to separate them cleanly into two independent source-to-target paths.
            
            result.path1.reserve(p1_edges.size());
            for (const auto& it : p1_edges) {
                result.path1.push_back(*it);
            }
            result.path2.reserve(p2_edges.size());
            for (const auto& it : p2_edges) {
                result.path2.push_back(*it);
            }
            result.success = true;
            return result;
        };
    };

}