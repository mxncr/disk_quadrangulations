#include <stdlib.h>
#include <iostream>
#include <numeric>
#include <limits.h>


#include "third_party/cxxopts.hpp"
#include "src/logging.h"
#include "src/disk_quadrangulation.h"
#include "src/cpp_utils.hpp"

using namespace Logging;
using std::string;
using std::vector;
using namespace DiskQuadrangulation;

int main(int argc, char** argv) {
  std::string output_filename;
  int qmax = INT_MAX;
  int nbv_bdr = 0;
  int nbv_int = 0;
  int limit = INT_MAX;
  bool print = false;
  bool no_symmetry = false;

  try {
    cxxopts::Options options("dquad", "dquad --- Compute topological quadrangulations of the disk");
    options
      .positional_help("<output.dquads>")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
      ("o,output", "Quadrangulations", cxxopts::value<std::string>()->default_value("output.dquads"),"NEW_FILE")
      ("q,max-quads", "Maximum number of quads in the results", cxxopts::value<int>(qmax),"INT")
      ("b,nb-boundary-vertices", "Number of vertices on the disk", cxxopts::value<int>(nbv_bdr),"INT")
      ("i,nb-interor-vertices", "Number of vertices strictly inside the disk", cxxopts::value<int>(nbv_int),"INT")
      ("l,limit-quadrangulations", "Stop the search after l quadrangulations are found", cxxopts::value<int>(limit),"INT")
      ("n,no-symmetry", "Do not use symmetries to reduce search", cxxopts::value<bool>(no_symmetry),"BOOL")
      ("p,print", "Print all the quadrangulations", cxxopts::value<bool>(print),"BOOL")
      ("h,help", "Print help")
      ;

    options.parse_positional({"output"});

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help({""}) << std::endl;
      exit(0);
    }
    output_filename = result["output"].as<std::string>();
  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }


  if (nbv_bdr > 0) {
    size_t vmax = nbv_int + nbv_bdr + 0;
    vector<id> disk(nbv_bdr);
    std::iota(disk.begin(), disk.end(), 0);
    info("input disk has {} vertices", nbv_bdr);
    info("computing the disk quadrangulations ...");
    info("exluding symmetry: {}", !no_symmetry);

    DiskQuadrangulationOpt opt;
    opt.nb_qdrl_max = limit;
    opt.nb_vert_max = nbv_int + nbv_bdr;
    opt.nb_quad_max = qmax;
    opt.symmetry = !no_symmetry;

    info("disk quadrangulations containing {} bdr. and {} interior vertices, with max {} quads and max {} vertices", nbv_bdr, nbv_int, qmax, vmax);
    vector<QuadMesh> meshes;
    find_quad_meshes(disk, opt, meshes);

    size_t maxq = 0;
    for (size_t i = 0; i < meshes.size(); ++i) {
      maxq = std::max(maxq,meshes[i].size());
    }

    info(" found {} quadrangulations, with {} quads at most", meshes.size(), maxq);
    if (print) {
      for (size_t i = 0; i < meshes.size(); ++i) {
        info("  {}: {}", i, meshes[i]);
      }
    }
    if (meshes.size() >= limit) {
      warn("reached search limit number of quadrangulations ({})", limit);
    }
  } else if (qmax != INT_MAX) {
    DiskQuadrangulationOpt opt;
    opt.nb_qdrl_max = limit;
    opt.nb_quad_max = qmax;
    opt.nb_vert_max = (id)(std::min(4+3*size_t(opt.nb_quad_max),(size_t)ID_MAX));
    opt.symmetry = !no_symmetry;

    info("generate disk quadrangulations with max {} quads", qmax);
    vector<QuadMesh> meshes;
    generate_quad_meshes(opt, meshes);
    size_t maxq = 0;
    size_t nvm = 0;
    for (size_t i = 0; i < meshes.size(); ++i) {
      for (size_t j = 0; j < meshes[i].size(); ++j) {
        for (size_t k = 0; k < meshes[i][j].size(); ++k) {
          nvm = std::max(nvm, (size_t) meshes[i][j][k]);
        }
      }
      maxq = std::max(maxq,meshes[i].size());
    }
    info(" found {} quadrangulations, with {} quads and {} vertices at most", meshes.size(), maxq, nvm+1);
    if (print) {
      for (size_t i = 0; i < meshes.size(); ++i) {
        size_t nvmi = 0;
        for (size_t j = 0; j < meshes[i].size(); ++j) {
          for (size_t k = 0; k < meshes[i][j].size(); ++k) {
            nvmi = std::max(nvmi, (size_t) meshes[i][j][k]);
          }
        }
        info("  {}: {}, nbv = {}", i, meshes[i], nvmi+1);
      }
    }
    if (meshes.size() >= limit) {
      warn("reached search limit number of quadrangulations ({})", limit);
    }

  }

  return 0;
}
