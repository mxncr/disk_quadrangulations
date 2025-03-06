#include <stdlib.h>
#include <iostream>
#include <numeric>
#include <optional>
#include <limits.h>

#include "src/logging.h"
#include "src/disk_quadrangulation.h"
#include "src/cpp_utils.hpp"

#include "third_party/cxxopts.hpp"

using namespace DiskQuadrangulation;
using namespace Logging;
using std::string;
using std::vector;


int main(int argc, char** argv) {
  DiskQuadrangulationOpt opt;

  try {
    cxxopts::Options options("enumerate", "enumerate --- Compute topological quadrangulations of the disk");
    options
      .positional_help("<output.dqs>")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
      ("o,output", "Quadrangulations", cxxopts::value<std::string>()->default_value(opt.output_file),"NEW_FILE")
      ("Q,nb-quad-max", "Maximum number of quads in the meshes", cxxopts::value<id>(opt.nb_quad_max),"ID")
      ("V,nb-vert-max", "Maximum number of vertices in the meshes", cxxopts::value<id>(opt.nb_vert_max),"ID")
      ("I,nb-int-vert-max", "Maximum number of interior vertices in the meshes", cxxopts::value<id>(opt.nb_int_vert_max),"ID")
      ("B,nb-bdr-vert-max", "Maximum number of boundary vertices in the meshes", cxxopts::value<id>(opt.nb_bdr_vert_max),"ID")
      ("val-int-max", "Maximum quad-valence of interior vertices", cxxopts::value<id>(opt.val_int_max),"ID")
      ("val-bdr-max", "Maximum quad-valence of boundary vertices", cxxopts::value<id>(opt.val_bdr_max),"ID")
      ("n,nb-qdrl-max", "Stop after N quadrangulations", cxxopts::value<size_t>(opt.nb_qdrl_max),"SIZE")
      ("v,verbose", "Print information", cxxopts::value<bool>(opt.verbose),"BOOL")
      ("viz", "Launch gmsh after", cxxopts::value<bool>(opt.visualization),"BOOL")
      ("cpp", "Print quad mesh arrays", cxxopts::value<bool>(opt.cpp),"BOOL")
      ("manual", "Custom selection of quadrangulations (for quadqs)", cxxopts::value<bool>(opt.manual),"BOOL")
      ("d,debug", "Print debug information", cxxopts::value<bool>(opt.debug),"BOOL")
      ("h,help", "Print help")
      ;

    options.parse_positional({"output"});

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help({""}) << std::endl;
      exit(0);
    }
    opt.output_file = result["output"].as<std::string>();
  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(1);
  }

  if (opt.debug) opt.verbose = true;
  Logging::LOG_ENABLE_TRACE = (opt.verbose && opt.debug);

  enumerate_disk_quadrangulations(opt);

  return 0;
}
