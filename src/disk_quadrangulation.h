#pragma once

#include <array>
#include <vector>
#include <string>
#include <limits>
#include <stdint.h>
#include <limits.h>

namespace DiskQuadrangulation {

  using id       = uint8_t;
  using id2      = std::array<id, 2>;
  using id4      = std::array<id, 4>;
  using DiskMesh = std::vector<id2>;
  using QuadMesh = std::vector<id4>;

  constexpr id NO_ID = std::numeric_limits<id>::max();
  constexpr id ID_MAX = std::numeric_limits<id>::max();

  struct DiskQuadrangulationOpt {
    id            nb_quad_max      = ID_MAX;
    id            nb_vert_max      = ID_MAX;
    id            nb_int_vert_max  = ID_MAX;
    id            nb_bdr_vert_max  = ID_MAX;
    id            val_int_max      = ID_MAX;
    id            val_bdr_max      = ID_MAX;
    size_t        nb_qdrl_max      = SIZE_MAX;
    bool          symmetry         = true;
    std::string   output_file      = "/tmp/output.dqs";
    bool          verbose          = true;
    bool          debug            = false;
    bool          visualization    = false;
    bool          no_nauty_check   = false;
    bool          cpp              = false;
    bool          manual           = false;
  };

  /* Find the quad meshes that match the bdr. disk */
  bool find_quad_meshes(
      const std::vector<id>& disk,
      const DiskQuadrangulationOpt& opt,
      std::vector<QuadMesh>& quadMeshes);

  /* Generate all quad meshes containing up to opt.nb_quad_max quads */
  bool generate_quad_meshes(
      const DiskQuadrangulationOpt& opt,
      std::vector<QuadMesh>& quadMeshes);

  bool enumerate_disk_quadrangulations(const DiskQuadrangulationOpt& opt);

}
