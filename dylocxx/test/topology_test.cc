

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



TEST_F(TopologyTest, DistanceMetric) {
  dyloc::init(&TESTENV.argc, &TESTENV.argv);

  auto & topo = dyloc::team_topology();

  const auto & dom_source = topo.domains().at(".0.0.0.0.0.0");
  const auto & dom_target = topo.domains().at(".0.0.0.0.0.0.0");

  auto dist_fn = [&](const auto & a,
                     const auto & b){
    return (7 - topo.ancestor({a.domain_tag, b.domain_tag}).level) *10; };

  std::string metric_name = "user_metric";
  topo.add_distance_metric(metric_name, dist_fn);

  DYLOC_LOG_DEBUG("TopologyTest.DistanceMetric", "Metrics:");
  auto dist_metrics = topo.list_distance_metrics();
  DYLOC_LOG_DEBUG("TopologyTest.DistanceMetric", dist_metrics);

  ASSERT_EQ(10, topo.calculate_distance(dom_source, dom_target, metric_name));
  
  std::for_each(topo.domains().begin(),
	        topo.domains().end(),
		[&](const auto & target){
    std::cout << "Distance ("<< metric_name <<") from " << dom_source.domain_tag << " to " << target.first << " >> " << std::to_string(topo.calculate_distance(dom_source,target.second, metric_name)) << std::endl;
                 });

  graphviz_out(topo.graph(), "testgraph.dot");

  dyloc::finalize();
}

TEST_F(TopologyTest, AddandRemoveConnection) {
  dyloc::init(&TESTENV.argc, &TESTENV.argv);
  auto & topo = dyloc::team_topology();
  DYLOC_LOG_DEBUG("connect()","start testing");
  const std::string source_tag = ".0";
  const std::string target_tag = ".0.0.0.0";

  DYLOC_LOG_DEBUG("attempting to connect",source_tag, target_tag);
  //same source/target but different types of edges
  auto edge1 = topo.connect(source_tag, target_tag, topology::edge_properties{topology::edge_type::contains, 42});
  auto edge2 = topo.connect(source_tag, target_tag, topology::edge_properties{topology::edge_type::logical, 43});
  ASSERT_EQ(1, edge1.second);
  ASSERT_EQ(1, edge2.second);
  auto edge3 = topo.connect(source_tag, target_tag, topology::edge_properties{topology::edge_type::contains, 42});
  ASSERT_NE(0, edge3.second);// 2 same connections? TODO!?
  DYLOC_LOG_DEBUG("connected sucessfully","Assertion completed");
  graphviz_out(topo.graph(), "testgraph1.dot");
  //-> disconnect
  topo.disconnect(edge1.first);
  topo.disconnect(edge2.first);
  topo.disconnect(edge3.first);

//ASSERT_EQ(topo,topo_cpy); would have to implement operator= (check all members for equality)
  graphviz_out(topo.graph(), "testgraph2.dot");

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

