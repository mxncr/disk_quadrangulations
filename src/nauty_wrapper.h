#pragma once

#include <vector>
#include <array>
#include <string>

namespace Nauty {

  using id = uint8_t;
  using id2 = std::array<uint8_t,2>;

  bool graph_canonical_labelling(int N, const std::vector<id2>& edges, std::vector<id2>& canonical_edges);

}
