#ifndef DYLOCXX__RUNTIME_H__INCLUDED
#define DYLOCXX__RUNTIME_H__INCLUDED

#include <dylocxx/host_topology.h>
#include <dylocxx/unit_mapping.h>
#include <dylocxx/locality_domain.h>
#include <dylocxx/domain_graph.h>

#include <dyloc/common/types.h>

#include <dash/dart/if/dart_types.h>

#include <vector>


namespace dyloc {

class runtime {
  std::unordered_map<dart_team_t, host_topology>   _host_topologies;
  std::unordered_map<dart_team_t, unit_mapping>    _unit_mappings;
  std::unordered_map<dart_team_t, locality_domain> _locality_domains;

  std::unordered_map<dart_team_t, domain_graph>    _domain_graphs;

 public:
  void initialize();
  void finalize();

  void initialize_locality(dart_team_t team);
  void finalize_locality(dart_team_t team);

  const dyloc_unit_locality_t & unit_locality(
          dart_team_t t,
          dart_team_unit_t u) {
    return _unit_mappings.at(t)[u];
  }

  const dyloc_unit_locality_t & unit_locality(
          dart_global_unit_t u) {
    // Unit id in team ALL is identical to global unit id:
    return _unit_mappings.at(DART_TEAM_ALL)[u.id];
  }

  domain_graph & locality_graph(
    dart_team_t t) {
    return _domain_graphs.at(t);
  }
};

} // namespace dyloc

#endif // DYLOCXX__RUNTIME_H__INCLUDED
