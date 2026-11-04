#pragma once
#include <functional>
#include <vector>
#include <ranges>
#include <algorithm>
#include <optional>

#include <type_traits>
#include <limits>
#include <graph/graph.hpp>

namespace amet{

    template<typename G>
    class vf2 {
        public:
        vf2() = default;
        vf2(const vf2& orig) = delete;
        virtual ~vf2() = default;

        void operator()(std::function<void()> apply){

        }

    };

    template <typename G>
    auto prepare_vf2_order(const G& g) {
        using vertex_id_t = graph::vertex_id_t<G>;
        std::vector<vertex_id_t> vertices;
        vertices.reserve(num_vertices(g));
        
        // Populate indices
        for(vertex_id_t i = 0; i < num_vertices(g); ++i) {
            vertices.push_back(i);
        }

        // Sort by a composite heuristic: Degree combined with coordination layout
        std::ranges::sort(vertices, [&g](vertex_id_t a, vertex_id_t b) {
            auto deg_a = graph::degree(g, a);
            auto deg_b = graph::degree(g, b);
            if (deg_a != deg_b) return deg_a > deg_b;
            
            // Secondary sort: sum of neighbor degrees (2-look-ahead density)
            auto sum_neighbors = [](const auto& graph, vertex_id_t v) {
                size_t count = 0;
                for (auto&& edge : graph::adj_list::edges(graph, v)) {
                    count += graph::degree(graph, edge.target_id());
                }
                return count;
            };
            return sum_neighbors(g, a) > sum_neighbors(g, b);
        });

        return vertices;
    }

    template <typename G>
    constexpr auto invalid_vertex_id() -> graph::vertex_id_t<G>{
        using vertex_id_t = graph::vertex_id_t<G>;

        if constexpr (std::is_integral_v<vertex_id_t>) {
            // Rule 1: Primitive integers (standard in dynamic_graph configs)
            return std::numeric_limits<vertex_id_t>::max();
        } 
        else if constexpr (std::is_default_constructible_v<vertex_id_t>) {
            // Rule 2: Custom structural descriptors / proxies
            // If graph-v3 transitions to a non-primitive layout type, its standard
            // default state typically acts as the internal library sentinel.
            return vertex_id_t{}; 
        }
        else {
            // Rule 3: Compile-time safeguard
            static_assert(sizeof(G) == 0, 
                "amet error: graph::vertex_id_t is a complex type with no visible default constructor.");
        }
    };

template <typename PatternG, typename TargetG>
struct vf2_state {

    const PatternG& pattern;
    const TargetG&  target;

    // Bijection maps: maps vertex from one graph to the other (or invalid_vertex_id)
    std::vector<graph::vertex_id_t<TargetG>>  core_p; 
    std::vector<graph::vertex_id_t<PatternG>> core_t;

    // Vectors to track the boundary/terminal sets for structural pruning
    std::vector<size_t> in_p;  // Depth at which vertex entered pattern boundary
    std::vector<size_t> in_t;  // Depth at which vertex entered target boundary
    
    size_t depth = 0;

    vf2_state(const PatternG& p, const TargetG& t) 
        : pattern(p), target(t),
          core_p(num_vertices(p), amet::invalid_vertex_id<PatternG>()),
          core_t(num_vertices(t), amet::invalid_vertex_id<TargetG>()),
          in_p(num_vertices(p), 0),
          in_t(num_vertices(t), 0) {}

};

template <typename PatternG, typename TargetG>
class VF2PlusPlus {
public:
    VF2PlusPlus() = default;
    virtual ~VF2PlusPlus() = default;

    // Core search loop using customizable projection predicates for chemistry logic
    template <typename VertexMatchPredicate, typename EdgeMatchPredicate, typename MatchCallback>
    void match(const PatternG& pattern, const TargetG& target, 
               VertexMatchPredicate v_match, EdgeMatchPredicate e_match, MatchCallback on_match) 
    {
        vf2_state state(pattern, target);
        std::vector<std::pair<graph::vertex_id_t<PatternG>, graph::vertex_id_t<TargetG>>> current_bijection;
        
        recursive_match(state, v_match, e_match, on_match, current_bijection);
    }

private:
    template <typename VertexMatchPredicate, typename EdgeMatchPredicate, typename MatchCallback, typename BijectionT>
    void recursive_match(vf2_state<PatternG, TargetG>& state, 
                         VertexMatchPredicate& v_match, EdgeMatchPredicate& e_match,
                         MatchCallback& on_match, BijectionT& current_bijection) 
    {
        if (state.depth == num_vertices(state.pattern)) {
            // Found a valid subgraph isomorphism! Record the bijection map.
            on_match(current_bijection);
            return;
        }

        // Generate candidate pairs using the VF2++ look-ahead ordering heuristic
        auto candidates = generate_candidates(state);

        for (const auto& [p_candidate, t_candidate] : candidates) {
            if (is_feasible(state, p_candidate, t_candidate, v_match, e_match)) {
                
                // Push state
                push_pair(state, p_candidate, t_candidate);
                current_bijection.emplace_back(p_candidate, t_candidate);

                recursive_match(state, v_match, e_match, on_match, current_bijection);

                // Pop state (Backtrack)
                current_bijection.pop_back();
                pop_pair(state, p_candidate, t_candidate);
            }
        }
    }

    template <typename VertexMatchPredicate, typename EdgeMatchPredicate>
    bool is_feasible(const vf2_state<PatternG, TargetG>& state, 
                    graph::vertex_id_t<PatternG> p_v, graph::vertex_id_t<TargetG> t_v,
                    VertexMatchPredicate& v_match, EdgeMatchPredicate& e_match) 
    {
        using p_vertex_t = graph::vertex_id_t<PatternG>;
        using t_vertex_t = graph::vertex_id_t<TargetG>;

        // 1. Semantic Check (e.g., node properties/labels)
        if (!v_match(state.pattern, p_v, state.target, t_v)) {
            return false;
        }

        // Structural Look-Ahead Counters
        size_t p_back_count = 0, t_back_count = 0; // Edges to Core
        size_t p_bnd_count  = 0, t_bnd_count  = 0; // Edges to Boundary (Terminal set)
        size_t p_rem_count  = 0, t_rem_count  = 0; // Edges to Remainder (Unvisited)

        // --- Examine Pattern Node Neighbors ---
        for (auto&& p_edge : graph::adj_list::edges(state.pattern, p_v)) {
            p_vertex_t p_neighbor = graph::target_id(state.pattern, p_edge);
            t_vertex_t t_neighbor = state.core_p[p_neighbor];

            if (t_neighbor != amet::invalid_vertex_id<TargetG>()) {
                // Rule: If neighbor is matched, the target node MUST have a matching edge to t_neighbor
                auto t_edges = graph::adj_list::edges(state.target, t_v);
                auto it = std::ranges::find_if(t_edges, [&state, t_neighbor](const auto& e) {
                    return graph::target_id(state.target, e) == t_neighbor;
                });

                if (it == std::ranges::end(t_edges) || !e_match(state.pattern, p_edge, state.target, *it)) {
                    return false; 
                }
                p_back_count++;
            } else if (state.in_p[p_neighbor] > 0) {
                p_bnd_count++;
            } else {
                p_rem_count++;
            }
        }

        // --- Examine Target Node Neighbors ---
        for (auto&& t_edge : graph::adj_list::edges(state.target, t_v)) {
            t_vertex_t t_neighbor = graph::target_id(state.target, t_edge);
            p_vertex_t p_neighbor = state.core_t[t_neighbor];

            if (p_neighbor != amet::invalid_vertex_id<PatternG>()) {
                t_back_count++;
            } else if (state.in_t[t_neighbor] > 0) {
                t_bnd_count++;
            } else {
                t_rem_count++;
            }
        }

        // 2. Look-Ahead Pruning Invariants (Crucial for exact Subgraph Isomorphism counts)
        // Target vertex must have at least as many connections to each structural zone as the pattern vertex
        if (p_bnd_count > t_bnd_count) return false; // 1-Look-Ahead Failure
        if (p_rem_count > t_rem_count) return false; // 2-Look-Ahead Failure

        return true;
    }

    // VF2++ optimization: local boundary tracking allocations go here...
    // 2. Commit a matched pair to the state matrices and expand boundaries
    void push_pair(vf2_state<PatternG, TargetG>& state, graph::vertex_id_t<PatternG> p_v, graph::vertex_id_t<TargetG> t_v) {
        state.depth++;
        state.core_p[p_v] = t_v;
        state.core_t[t_v] = p_v;

        // Update boundary depths for neighbors of p_v
        for (auto&& edge : graph::adj_list::edges(state.pattern, p_v)) {
            auto neighbor = graph::target_id(state.pattern, edge);
            if (state.in_p[neighbor] == 0) state.in_p[neighbor] = state.depth;
        }

        // Update boundary depths for neighbors of t_v
        for (auto&& edge : graph::adj_list::edges(state.target, t_v)) {
            auto neighbor =  graph::target_id(state.target, edge);
            if (state.in_t[neighbor] == 0) state.in_t[neighbor] = state.depth;
        }
    }

    // 3. Roll back a state during backtracking
    void pop_pair(vf2_state<PatternG, TargetG>& state, graph::vertex_id_t<PatternG> p_v, graph::vertex_id_t<TargetG> t_v) {
        state.core_p[p_v] = amet::invalid_vertex_id<TargetG>();
        state.core_t[t_v] = amet::invalid_vertex_id<PatternG>();

        for (graph::vertex_id_t<PatternG> i = 0; i < num_vertices(state.pattern); ++i) {
            if (state.in_p[i] == state.depth) state.in_p[i] = 0;
        }
        for (graph::vertex_id_t<TargetG> i = 0; i < num_vertices(state.target); ++i) {
            if (state.in_t[i] == state.depth) state.in_t[i] = 0;
        }
        state.depth--;
    }

    auto generate_candidates(const vf2_state<PatternG, TargetG>& state)-> std::vector<std::pair<graph::vertex_id_t<PatternG>, graph::vertex_id_t<TargetG>>> {
        using p_vertex_t = graph::vertex_id_t<PatternG>;
        using t_vertex_t = graph::vertex_id_t<TargetG>;
        std::vector<std::pair<p_vertex_t, t_vertex_t>> candidates;

        // Find the first unmatched pattern vertex that belongs to the pattern boundary
        std::optional<p_vertex_t> min_p;
        for (p_vertex_t i = 0; i < num_vertices(state.pattern); ++i) {
            if (state.core_p[i] == amet::invalid_vertex_id<TargetG>() && state.in_p[i] > 0) {
                min_p = i;
                break;
            }
        }

        // If the boundary set is empty, pick the first completely unmapped pattern vertex
        if (!min_p.has_value()) {
            for (p_vertex_t i = 0; i < num_vertices(state.pattern); ++i) {
                if (state.core_p[i] == amet::invalid_vertex_id<TargetG>()) {
                    min_p = i;
                    break;
                }
            }
        }

        if (!min_p.has_value()) return candidates; // Safety catch

        p_vertex_t p_v = *min_p;

        // Pair it with unmatched target vertices on the target boundary (or any unmatched target if boundary empty)
        bool target_boundary_has_nodes = std::ranges::any_of(state.in_t, [](size_t depth) { return depth > 0; });

        for (t_vertex_t t_v = 0; t_v < num_vertices(state.target); ++t_v) {
            if (state.core_t[t_v] == amet::invalid_vertex_id<PatternG>()) {
                if (target_boundary_has_nodes && state.in_t[t_v] == 0) {
                    continue; // Skip target vertices not on the boundary if a boundary exists
                }
                candidates.emplace_back(p_v, t_v);
            }
        }

        return candidates;
    }
};

}