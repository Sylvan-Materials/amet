#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <string>
#include <vector>
#include <array>
#include <functional>

#include <graph/graph.hpp>
#include "graph/container/dynamic_graph.hpp"
#include <graph/container/traits/vov_graph_traits.hpp>
#include "graph/adaptors/filtered_graph.hpp" // Adaptor to prune branches
#include "graph/views/incidence.hpp"
#include "graph/views/vertexlist.hpp"
#include "graph/views/edgelist.hpp"

#include <amet/automorphism.h>
#include <amet/Suurballe.h>
#include <amet/VF2.h>
#include <amet/connected_components.h>
#include <amet/biconnected_components.h>

// Explicitly modeling a minimal graph-v3 compliant structural array
struct MockGraph {
    using vertex_descriptor = size_t;
    using edge_descriptor = size_t;

    std::vector<size_t> vertices{0, 1, 2, 3};
    
    auto begin() const { return vertices.begin(); }
    auto end() const { return vertices.end(); }
};

// High-precision representation of the Petersen Graph
struct PetersenGraph {
    using vertex_descriptor = size_t;
    
    // 10 Vertices, 3-regular connectivity matrix
    static constexpr size_t num_vertices = 10;
    static constexpr std::array<std::array<size_t, 3>, 10> adjacency_list{{
        {1, 4, 5}, {0, 2, 6}, {1, 3, 7}, {2, 4, 8}, {3, 0, 9}, // Outer pentagon
        {0, 7, 8}, {1, 8, 9}, {2, 5, 9}, {3, 5, 6}, {4, 6, 7}  // Inner star
    }};
};

// Benzene molecule topological wireframe (excluding explicit coordinates)
struct BenzeneGraph {
    using vertex_descriptor = size_t;
    static constexpr size_t num_vertices = 12;
    
    // Indices 0-5 are Carbon, 6-11 are matching bonded Hydrogens
    static constexpr std::array<std::string_view, 12> element_colors{
        "C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"
    };
};

// Define a tiny, generic tuple-like record structure directly in the test file 
// to satisfy graph-v3's structured binding vertex requirements natively
struct abstract_vertex_record {
    unsigned int id;
    int value;

    template <std::size_t I>
    decltype(auto) get() const {
        if constexpr (I == 0) return id;
        if constexpr (I == 1) return value;
    }
};

namespace std {
    template <> struct tuple_size<abstract_vertex_record> : integral_constant<size_t, 2> {};
    template <size_t I> struct tuple_element<I, abstract_vertex_record> {
        using type = conditional_t<I == 0, unsigned int, int>;
    };
}

// Simple, primitive edge record structure for standard developer tests
struct abstract_edge_record {
    unsigned int source_id;
    unsigned int target_id;
    int value;
};

using G = graph::container::dynamic_adjacency_graph<graph::container::vov_graph_traits<std::string, std::string>>;


TEST_SUITE("algorithms"){
    
TEST_CASE("Validating amet Graph Type Constraints") {
    MockGraph g;
    
    // Validate that our C++23 viewable_incidence_graph concept holds true
    REQUIRE(amet::viewable_incidence_graph<MockGraph>);
}

TEST_CASE("Validating Partition Refiner Invariants") {
    MockGraph g;
    
    // Ensure that our zero-allocation projection layer compiles cleanly
    auto refiner = amet::partition_refiner(g, [](size_t v) { return v % 2; });
    auto orbits = refiner.compute_orbits();
    
    CHECK(orbits.size() == 0); // Update when your main splitting logic drops in!
}

TEST_CASE("Validate Suurballe shortest path") {
    G diamondGraph;
    std::vector<graph::copyable_vertex_t<graph::vertex_id_t<G>, std::string>> bases;
    std::vector<graph::copyable_edge_t<graph::vertex_id_t<G>, std::string>> paths;
    bases.push_back({bases.size(), "First"});
    bases.push_back({bases.size(), "Second"});
    bases.push_back({bases.size(), "Third"});
    bases.push_back({bases.size(), "Home"});
    paths.push_back({0, 1, "First-to-Second"});
    paths.push_back({0, 2, "First-to-Third"});
    paths.push_back({1, 2, "Second-to-Third"});
    paths.push_back({2, 3, "Third-to-Home"});
    paths.push_back({3, 0, "Home-to-First"});

    diamondGraph.reserve_vertices(bases.size());
    diamondGraph.load_vertices(bases, std::identity{});

    diamondGraph.reserve_edges(paths.size());
    diamondGraph.load_edges(paths, std::identity{});

    auto uniform_weight = [](const auto& gr, const auto& uv) -> double { 
        return 1.0; 
    };

    std::cout<<"suurballe "<<std::endl;
    amet::Suurballe<G> suurballe;
    graph::vertex_id_t<G> homeVertex(3);
    auto secondVertex = graph::vertex_id_t<G>(1);
    amet::suurballe_result<G> result=suurballe(diamondGraph, secondVertex, homeVertex, uniform_weight);
    std::cout<<"result.success "<<result.success<<std::endl;
    REQUIRE(result.success);
    CHECK(result.path1.size() == 2);
    CHECK(result.path2.size() == 2);
}

TEST_CASE("VF2++ Subgraph Isomorphism: Diamond in a Wheel Graph") {
    // We define a simple graph type using graph-v3 structures
    // Vertices store a string label, Edges are string as well, Allocator is standard
    using UndirectedGraph = graph::container::dynamic_adjacency_graph<graph::container::vov_graph_traits<std::string, std::string, void, unsigned int, true>>;

    // =========================================================================
    // 1. Construct the Target UndirectedGraph: Wheel UndirectedGraph W5 (5 vertices, 8 edges)
    // =========================================================================
    UndirectedGraph wheel;
    std::vector<graph::copyable_vertex_t<graph::vertex_id_t<UndirectedGraph>, std::string>> wheel_nodes;
    std::vector<graph::copyable_edge_t<graph::vertex_id_t<UndirectedGraph>, std::string>> wheel_edges;
    // Hub is index 0, Rim nodes are 1, 2, 3, 4
    wheel_nodes.push_back({0, "Hub"});
    wheel_nodes.push_back({1, "Rim 1"});
    wheel_nodes.push_back({2, "Rim 2"});
    wheel_nodes.push_back({3, "Rim 3"});
    wheel_nodes.push_back({4, "Rim 4"});

    // Spokes connecting hub to rim
    wheel_edges.push_back({0, 1, "Hub-to-Rim1"});
    wheel_edges.push_back({0, 2, "Hub-to-Rim2"});
    wheel_edges.push_back({0, 3, "Hub-to-Rim3"});
    wheel_edges.push_back({0, 4, "Hub-to-Rim4"});

    // Outer rim cycle

    // Rim 1
    wheel_edges.push_back({1, 2, "Rim1-to-Rim2"});

    // Rim 2
    wheel_edges.push_back({2, 3, "Rim2-to-Rim3"});

    // Rim 3
    wheel_edges.push_back({3, 4, "Rim3-to-Rim4"});

    // Rim 4
    wheel_edges.push_back({4, 1, "Rim4-to-Rim1"});

    wheel.reserve_vertices(wheel_nodes.size());
    wheel.load_vertices(wheel_nodes, std::identity{});

    wheel.reserve_edges(wheel_edges.size());
    wheel.load_edges(wheel_edges, std::identity{});

    // =========================================================================
    // 2. Construct the Pattern UndirectedGraph: Diamond UndirectedGraph (4 vertices, 5 edges)
    // =========================================================================
    UndirectedGraph diamond;
    std::vector<graph::copyable_vertex_t<graph::vertex_id_t<UndirectedGraph>, std::string>> diamond_nodes;
    std::vector<graph::copyable_edge_t<graph::vertex_id_t<UndirectedGraph>, std::string>> diamond_edges;
    diamond_nodes.push_back({diamond_nodes.size(), "Top"});        // ID: 0
    diamond_nodes.push_back({diamond_nodes.size(), "Left"});       // ID: 1
    diamond_nodes.push_back({diamond_nodes.size(), "Right"});      // ID: 2
    diamond_nodes.push_back({diamond_nodes.size(), "Bottom"});     // ID: 3

    // Source 1
    diamond_edges.push_back({0, 1, "Top-to-Left"});
    diamond_edges.push_back({0, 2, "Top-to-Right"});
    diamond_edges.push_back({0, 3, "Top-to-Bottom"});
    
    // Source 1
    diamond_edges.push_back({1, 3, "Left-to-Bottom"});

    // Source 2
    diamond_edges.push_back({2, 3, "Right-to-Bottom"});

    diamond.reserve_vertices(diamond_nodes.size());
    diamond.load_vertices(diamond_nodes, std::identity{});

    diamond.reserve_edges(diamond_edges.size());
    diamond.load_edges(diamond_edges, std::identity{});

    // =========================================================================
    // 3. Execution and Assertions
    // =========================================================================
    amet::VF2PlusPlus<UndirectedGraph, UndirectedGraph> matcher;
    
    size_t match_count = 0;
    
    // Simple pass-through predicates for pure topological matching
    auto vertex_always_match = [](const auto&, auto, const auto&, auto) { return true; };
    auto edge_always_match   = [](const auto&, const auto&, const auto&, const auto&) { return true; };

    // Run the matching loop, passing a lambda to capture solution discoveries
    matcher.match(diamond, wheel, vertex_always_match, edge_always_match, 
                  [&match_count](const auto& bijection) {
                      match_count++;
                      // Optional: verify bijection map invariants
                      CHECK(bijection.size() == 4); 
                  });

    // A Diamond graph embeds into a Wheel-5 exactly 8 times due to symmetry
    CHECK(match_count == 16);
}

using TestCG = graph::container::dynamic_adjacency_graph<graph::container::vov_graph_traits<int, int, void, unsigned int, true>>;

TEST_CASE("amet::connected_components - Standard Disconnected Clusters Benchmark") {
    // Topology Layout: 
    // Cluster 0 (Triangle Loop):   Nodes 0-1-2-0
    // Cluster 1 (Star / Tree):     Nodes 3-4, 3-5, 3-6
    // Cluster 2 (Isolated Node):  Node 7
    // Total Vertices: 8 (0 to 7)

    TestCG test_graph;
    
    // 1. Populate Vertices via abstract, local primitive records
    std::vector<abstract_vertex_record> vertices;
    for (unsigned int i = 0; i < 8; ++i) {
        vertices.push_back(abstract_vertex_record{.id = i, .value = 0});
    }
    
    test_graph.load_vertices(vertices, [](const abstract_vertex_record& v) {
        return std::make_pair(v.id, v.value);
    });

    // 2. Populate Edges Bidirectionally
    std::vector<abstract_edge_record> edges;
    auto add_edge = [&](unsigned int u, unsigned int v) {
        edges.push_back(abstract_edge_record{.source_id = u, .target_id = v, .value = 1});
        edges.push_back(abstract_edge_record{.source_id = v, .target_id = u, .value = 1});
    };

    // --- Cluster 0: Triangle Loop ---
    add_edge(0, 1);
    add_edge(1, 2);
    add_edge(2, 0);

    // --- Cluster 1: Star Configuration ---
    add_edge(3, 4);
    add_edge(3, 5);
    add_edge(3, 6);

    // --- Cluster 2: Node 7 remains completely isolated (No edges added) ---

    test_graph.load_edges(edges, std::identity{});

    // 3. Execute your clean abstract amet components engine
    std::map<unsigned int, std::size_t> component_map;
    std::size_t total_components = amet::connected_components(test_graph, component_map);

    // 4. Assert Structural Invariants Rigorously
    // We expect exactly 3 separate isolated networks
    REQUIRE(total_components == 3);

    // Verify Cluster 0 (Triangle nodes share the same label)
    CHECK(component_map[0] == component_map[1]);
    CHECK(component_map[1] == component_map[2]);

    // Verify Cluster 1 (Star nodes share the same label)
    CHECK(component_map[3] == component_map[4]);
    CHECK(component_map[4] == component_map[5]);
    CHECK(component_map[5] == component_map[6]);

    // Verify absolute boundaries: The three clusters must possess unique IDs
    CHECK(component_map[0] != component_map[3]);
    CHECK(component_map[3] != component_map[7]);
    CHECK(component_map[0] != component_map[7]);
}

TEST_CASE("amet::biconnected_components - Standard Developer Benchmark Topology") {
    TestCG test_graph;
    
    // 1. Populate Vertices via abstract, local primitive records
    std::vector<abstract_vertex_record> vertices;
    for (unsigned int i = 0; i < 8; ++i) {
        vertices.push_back(abstract_vertex_record{.id = i, .value = 0});
    }
    
    // 3. THE FIX: Pass a projection lambda to explicitly feed the internal structured binding!
    test_graph.load_vertices(vertices, [](const abstract_vertex_record& v) {
        return std::make_pair(v.id, v.value);
    });

    // 4. Populate Edges Bidirectionally
    std::vector<abstract_edge_record> edges;
    auto add_edge = [&](unsigned int u, unsigned int v) {
        edges.push_back(abstract_edge_record{.source_id = u, .target_id = v, .value = 1});
        edges.push_back(abstract_edge_record{.source_id = v, .target_id = u, .value = 1});
    };

    // Block 1: Left Simple Cycle (0-1-3-2-0)
    add_edge(0, 1); add_edge(1, 3); add_edge(3, 2); add_edge(2, 0);

    // Central Cut-Edge (Bridge)
    add_edge(3, 5);

    // Block 2: Right Fused Cycle / Diamond Block (Nodes 3, 4, 5, 6)
    add_edge(5, 4); add_edge(4, 6); add_edge(6, 5); add_edge(3, 4);

    // Block 3: Isolated Leaf Spoke
    add_edge(5, 7);

    test_graph.load_edges(edges, std::identity{});
 
    // 3. Execute your clean amet engine
    std::vector<std::vector<unsigned int>> components;
    amet::biconnected_components(test_graph, components);

    // 4. Invariant Assertions
    REQUIRE(components.size() == 3);

    std::ranges::sort(components, [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });

    CHECK(components[0].size() == 4);
    CHECK(components[1].size() == 4);
    CHECK(components[2].size() == 2);
}

TEST_CASE("amet::biconnected_components - Fused and Isolated Molecular Ring System") {
    // Topology Layout: Naphthalene (fused 10-atom system) 
    // connected via a single bridge bond at Node 4 to a Phenyl Ring (6-atom system).
    // Total nodes: 16 (0 to 15).
    
    TestCG molecule_graph;
    
    // 1. Ingest Vertices using the clean pair-projection strategy
    std::vector<abstract_vertex_record> vertices;
    for (unsigned int i = 0; i < 16; ++i) {
        vertices.push_back(abstract_vertex_record{.id = i, .value = 0});
    }
    molecule_graph.load_vertices(vertices, [](const abstract_vertex_record& v) {
        return std::make_pair(v.id, v.value);
    });

    // 2. Ingest Bonds Bidirectionally
    std::vector<abstract_edge_record> edges;
    auto add_bond = [&](unsigned int u, unsigned int v) {
        edges.push_back(abstract_edge_record{.source_id = u, .target_id = v, .value = 1});
        edges.push_back(abstract_edge_record{.source_id = v, .target_id = u, .value = 1});
    };

    // --- Fused System (Naphthalene Block: Nodes 0 to 9) ---
    // Outer perimeter cycle
    add_bond(0, 1); add_bond(1, 2); add_bond(2, 3); add_bond(3, 4);
    add_bond(4, 5); add_bond(5, 6); add_bond(6, 7); add_bond(7, 8);
    add_bond(8, 9); add_bond(9, 0);
    // The central fused ring bridge bond sharing two interior carbons
    add_bond(1, 6); 

    // --- The Molecular Interfacial Bridge ---
    // Atom 4 is our structural Articulation Point! 
    add_bond(4, 10);

    // --- Isolated Pendant Ring (Phenyl Block: Nodes 10 to 15) ---
    add_bond(10, 11); add_bond(11, 12); add_bond(12, 13);
    add_bond(13, 14); add_bond(14, 15); add_bond(15, 10);

    molecule_graph.load_edges(edges, std::identity{});

    // 3. Execute the pure abstract amet engine
    std::vector<std::vector<unsigned int>> components;
    amet::biconnected_components(molecule_graph, components);

    // 4. Invariant Assertions
    // Unlike the dense diamond graph, these systems are separated by a true 
    // structural path bridge (4, 10) that is completely outside the cycles.
    // - Component 1: The Fused Naphthalene block (10 vertices)
    // - Component 2: The Isolated Phenyl ring (6 vertices)
    // - Component 3: The single connecting Bridge link {4, 10} (2 vertices)
    // Total Blocks = 3.
    REQUIRE(components.size() == 3);

    // Sort subsets by size descending for clean tracking
    std::ranges::sort(components, [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });

    CHECK(components[0].size() == 10); // Naphthalene block
    CHECK(components[1].size() == 6);  // Phenyl block
    CHECK(components[2].size() == 2);  // {4, 10} Interfacial Bridge

    // Helper lambda to check presence cleanly
    auto holds_vertex = [](const std::vector<unsigned int>& comp, unsigned int id) -> bool {
        return std::ranges::find(comp, id) != comp.end();
    };

    // Explicitly verify that Node 4 acts as the true articulation point
    // by ensuring it is present in both the Fused Block and the Bridge link.
    size_t node_4_occurrences = 0;
    for (const auto& comp : components) {
        if (holds_vertex(comp, 4)) {
            node_4_occurrences++;
        }
    }
    CHECK(node_4_occurrences == 2);
}

}
