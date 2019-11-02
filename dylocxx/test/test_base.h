#ifndef DYLOCXX__TEST__TEST_BASE_H__INCLUDED
#define DYLOCXX__TEST__TEST_BASE_H__INCLUDED

#include <gtest/gtest.h>

#include <dylocxx/topology.h>
#include <dylocxx/unit_locality.h>
#include <dylocxx/locality_domain.h>
#include <dylocxx/hwinfo.h>

#include <dylocxx/utility.h>
#include <dylocxx/adapter/dart.h>
#include <dylocxx/internal/logging.h>

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/depth_first_search.hpp>

//#if DYLOC_ENABLE_GRAPHVIZ
#include <boost/graph/graphviz.hpp>
//#endif

#include <string>
#include <iostream>
#include <iomanip>


namespace dyloc {
namespace test {

template <class DomainMap>
class locality_domain_dfs_output_visitor
: public boost::default_dfs_visitor {
  const DomainMap & _domains;
public:
  locality_domain_dfs_output_visitor(
    const DomainMap & domain_map)
  : _domains(domain_map) { }

  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, const Graph & g) const {
    const std::string            & ldom_tag = g[u].domain_tag;
    const auto &                   ldom_it  = _domains.find(ldom_tag);
    if (ldom_it == _domains.end()) {
      DYLOC_LOG_ERROR("locality_domain_dfs_output_visitor.discover_vertex",
                      "domain not found:", ldom_tag);
      return;
    }
    std::string ldom_st
      = ( g[u].state == dyloc::topology::vertex_state::hidden
        ? "hidden"
        : g[u].state == dyloc::topology::vertex_state::selected
        ? "selected"
        : g[u].state == dyloc::topology::vertex_state::moved
        ? "moved"
        : "" );
    const dyloc::locality_domain & ldom     = ldom_it->second;
    auto  ldom_maxthreads                   = ldom.host.size() > 0 &&
                                              ldom.team != DART_TEAM_NULL
                                              ? ldom.hwinfo().max_threads
                                              : -1;
    auto  ldom_shcache                      = ldom.host.size() > 0 &&
                                              ldom.team != DART_TEAM_NULL
                                              ? ldom.hwinfo().cache_sizes[2] /
                                                1024
                                              : -1;

    std::ostringstream ldom_sc_s;
    ldom_sc_s << std::string(ldom.level, ' ')
              << ldom.scope;

    std::cout << std::left  << std::setw(4)  << u               << " "
              << std::left  << std::setw(12) << ldom.host       << " "
              << std::left  << std::setw(9)  << ldom_st         << " "
              << std::left  << std::setw(14) << ldom_sc_s.str() << " "
              << std::left  << std::setw(4)  << ldom.level      << " g:"
              << std::right << std::setw(2)  << ldom.g_index    << " |"
              << std::left  << std::setw(20) << ldom_tag        << " |"
              << " nc:"  << std::right << std::setw(3) << ldom.num_cores  
              << " nt:"  << std::right << std::setw(3) << ldom_maxthreads
              << " l3d:" << std::right << std::setw(5) << ldom_shcache << "K"
              << " "
              << "units:"
              << dyloc::make_range(
                  ldom.unit_ids.begin(),
                  ldom.unit_ids.end())
              << std::endl;
  }
};

template <class Graph>
void graphviz_out(const Graph & graph, const std::string & filename) {
//#if DYLOC_ENABLE_GRAPHVIZ
  std::ofstream of(filename);
  write_graphviz(of, graph,
                 boost::make_label_writer(
                   boost::get(
                     &topology::vertex_properties::domain_tag,
                     graph)),
                 boost::make_label_writer(
                   boost::get(
                     &topology::edge_properties::distance,
                     graph)));
//#endif
}


class TestBase : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }
};

} // namespace test
} // namespace dyloc

#endif // DYLOCXX__TEST__TEST_BASE_H__INCLUDED
