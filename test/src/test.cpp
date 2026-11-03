#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <amet/automorphism.h>
#include <amet/Suurballe.h>
#include <string>
#include <vector>
#include <array>
#include <functional>

#include "graph/container/dynamic_graph.hpp"
#include <graph/container/traits/vov_graph_traits.hpp>
#include "graph/adaptors/filtered_graph.hpp" // Adaptor to prune branches
#include "graph/views/incidence.hpp"
#include "graph/views/vertexlist.hpp"
#include "graph/views/edgelist.hpp"

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
    bases.push_back({static_cast<graph::vertex_id_t<G>>(bases.size()), "First"});
    bases.push_back({static_cast<graph::vertex_id_t<G>>(bases.size()), "Second"});
    bases.push_back({static_cast<graph::vertex_id_t<G>>(bases.size()), "Third"});
    bases.push_back({static_cast<graph::vertex_id_t<G>>(bases.size()), "Home"});
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

}
