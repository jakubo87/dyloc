

#include <algorithm>
#include <memory>

#include "topology_test.h"
#include "test_globals.h"

#include <dylocxx/init.h>
#include <dylocxx/unit_locality.h>
#include <dylocxx/hwinfo.h>

#include <dylocxx/utility.h>
#include <dylocxx/adapter/dart.h>
#include <dylocxx/internal/logging.h>

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/depth_first_search.hpp>


namespace dyloc {
namespace test {


// Jakub
TEST_F(TopologyTest, DistanceMetric) {
  dyloc::init(&TESTENV.argc, &TESTENV.argv);

  auto & topo = dyloc::team_topology();

  //we want an iterator, dont we?
  //boost visitor?
  //
  const auto & dom_source_it = domain_vertices(".0.0");
  const auto & dom_target_it = domain_vertices(".0.1");
  //alternatively topo.domains().find(tag);

  DYLOC_LOG_DEBUG("TopologyTest.DistanceMetric", "Metrics:");

  auto dist_fn = [](const auto & a, const auto & b){ return 10; };
  std::string metric_name = "user_metric";
  topo.add_distance_metric(metric_name, dist_fn);

  auto dist_metrics = topo.list_distance_metrics();

  ASSERT_EQ(10, topo.calculate_distance(dom_source, dom_target, metric_name));

  DYLOC_LOG_DEBUG_VAR("TopologyTest.DistanceMetric", dist_metrics);

  dyloc::finalize();
}


TEST_F(TopologyTest, AddandRemoveConnection) {
  dyloc::init(&TESTENV.argc, &TESTENV.argv);
  auto & topo = dyloc::team_topology();

  const auto & source = domain_vertices[".0.0"];
  const auto & target = domain_vertices[".0.1"];

  //specifying category of edge
  using canon = topology::edge_type::contains;
  using log = topology::edge_type::logical;
  //same vertices but different types of edges
  auto edge1 = topo.connect(source, target, canon, 42);
  auto edge2 = topo.connect(source, target, log, 43);

  ASSERT_EQ(42, topo.get_edge_distance(edge1, canon));
  ASSERT_EQ(43, topo.get_edge_distance(edge2, log));

  //(hopefully) overwriting the edge
  auto edge2 = topo.connect(source, target, canon, 43);
  ASSERT_EQ(43, topo.get_edge_distance(edge1, canon));


  DYLOC_LOG_DEBUG_VAR("TopologyTest.AddandRemoveConnection");

  dyloc::finalize();
}

//TEST_F(TopologyTest, ExcludeDomains) {
//  dyloc::init(&TESTENV.argc, &TESTENV.argv);
//
//  auto & topo = dyloc::team_topology();
//  locality_domain_dfs_output_visitor<typename topology::domain_map>
//    vis(topo.domains());
//
//  DYLOC_LOG_DEBUG("TopologyTest.ExcludeDomains", "total domain hierarchy:");
//  if (dyloc::myid().id == 0) {
//    topo.depth_first_search(vis);
//  }
//  dart_barrier(DART_TEAM_ALL);
//
//  std::vector<std::string> excluded_domain_tags;
//  excluded_domain_tags.push_back(".0.0");
//
//  topo.exclude_domains(
//    excluded_domain_tags.begin(),
//    excluded_domain_tags.end());
//
//  if (dyloc::myid().id == 0) {
//    DYLOC_LOG_DEBUG("TopologyTest.ExcludeDomains",
//                    "domain hierarchy after excludes:");
//    topo.depth_first_search(vis);
//  }
//
//  ASSERT_EQ(0, topo.domains().count(".0.0"));
//
//  dyloc::finalize();
//}
//
//TEST_F(TopologyTest, DomainsAncestor) {
//  dyloc::init(&TESTENV.argc, &TESTENV.argv);
//
//  auto & topo = dyloc::team_topology();
//  locality_domain_dfs_output_visitor<typename topology::domain_map>
//    vis(topo.domains());
//
//  DYLOC_LOG_DEBUG("TopologyTest.DomainsAncestor", "total domain hierarchy:");
//  if (dyloc::myid().id == 0) {
//    topo.depth_first_search(vis);
//  }
//  dart_barrier(DART_TEAM_ALL);
//
//  int myid = dyloc::myid().id;
//
//  // Find ancestor domain from this vector of unit ids:
//  //
//  std::vector<dart_global_unit_t> unit_ids;
//  unit_ids.push_back(dyloc::myid());
//  for (int nd = 0; nd < std::min<int>(2, topo.domains().size()); ++nd) {
//    std::srand((nd + 1) * (myid + 1));
//    int random_unit_id = std::rand() % dyloc::num_units();
//    unit_ids.push_back(random_unit_id);
//  }
//
//  // Unit ids to domain tags:
//  //
//  std::vector<std::string> unit_domain_tags(unit_ids.size());
//  std::transform(unit_ids.begin(), unit_ids.end(),
//                 unit_domain_tags.begin(),
//                 [&](const dart_global_unit_t uid) -> std::string {
//                   auto unit_dom_tag = topo[uid].domain_tag;
//                   DYLOC_LOG_DEBUG_VAR("TopologyTest.DomainsAncestor",
//                                       unit_dom_tag);
//                   return unit_dom_tag;
//                 });
//
//  // Ancestor from unit domain tags:
//  //
//  const auto & ancestor = topo.ancestor(
//                            unit_domain_tags.begin(),
//                            unit_domain_tags.end());
//  for (const dart_global_unit_t & uid : unit_ids) {
//    ASSERT_NE(ancestor.unit_ids.end(),
//              std::find(ancestor.unit_ids.begin(),
//                        ancestor.unit_ids.end(),
//                        uid));
//  }
//  DYLOC_LOG_DEBUG_VAR("TopologyTest.DomainsAncestor", ancestor);
//
//  dyloc::finalize();
//}
//
//TEST_F(TopologyTest, ScopeDomains) {
//  dyloc::init(&TESTENV.argc, &TESTENV.argv);
//
//  auto & topo = dyloc::team_topology();
//  locality_domain_dfs_output_visitor<typename topology::domain_map>
//    vis(topo.domains());
//
//  DYLOC_LOG_DEBUG("TopologyTest.ScopeDomains", "total domain hierarchy:");
//  if (dyloc::myid().id == 0) {
//    topo.depth_first_search(vis);
//  }
//  dart_barrier(DART_TEAM_ALL);
//
//  // Tags of locality domains in CORE scope:
//  auto core_domain_tags = topo.scope_domain_tags(
//                            DYLOC_LOCALITY_SCOPE_CORE);
//  ASSERT_LT(0, core_domain_tags.size());
//  for (const auto & core_domain_tag : core_domain_tags) {
//    DYLOC_LOG_DEBUG_VAR("TopologyTest.ScopeDomains", core_domain_tag);
//    auto core_dom_it = topo.domains().find(core_domain_tag);
//    ASSERT_NE(topo.domains().end(), core_dom_it);
//    ASSERT_EQ(DYLOC_LOCALITY_SCOPE_CORE, core_dom_it->second.scope);
//  }
//
//  // Tags of locality domains in NUMA scope:
//  auto numa_domain_tags = topo.scope_domain_tags(
//                            DYLOC_LOCALITY_SCOPE_NUMA);
//  ASSERT_LT(0, numa_domain_tags.size());
//  for (const auto & numa_domain_tag : numa_domain_tags) {
//    DYLOC_LOG_DEBUG_VAR("TopologyTest.ScopeDomains", numa_domain_tag);
//    auto numa_dom_it = topo.domains().find(numa_domain_tag);
//    ASSERT_NE(topo.domains().end(), numa_dom_it);
//    ASSERT_EQ(DYLOC_LOCALITY_SCOPE_NUMA, numa_dom_it->second.scope);
//  }
//
//  dyloc::finalize();
//}



} // namespace dyloc
} // namespace test

