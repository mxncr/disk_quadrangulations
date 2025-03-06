#include "src/nauty_wrapper.h"
#include "src/logging.h" 
#include <iostream>
#include <fstream>

#include "third_party/nauty/nauty.h"

namespace Nauty {
  bool graph_canonical_labelling(int n, const std::vector<id2>& edges, std::vector<id2>& canonical_edges) {
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(int,map,map_sz);
    DYNALLSTAT(graph,g1,g1_sz);
    DYNALLSTAT(graph,cg1,cg1_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    options.getcanon = TRUE; // Select option for canonical labelling
    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC1(int,lab,lab_sz,n,"malloc");
    DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
    DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
    DYNALLOC1(int,map,map_sz,n,"malloc");
    DYNALLOC2(graph,g1,g1_sz,n,m,"malloc");
    DYNALLOC2(graph,cg1,cg1_sz,n,m,"malloc");

    // construct the nauty graph:
    EMPTYGRAPH(g1,m,n);
    for (size_t i = 0; i < edges.size(); ++i) {
      ADDONEEDGE(g1,edges[i][0],edges[i][1],m);
    }

    // get a canonical labeling of it:
    densenauty(g1,lab,ptn,orbits,&options,&stats,m,n,cg1);


    // extract the canonical edges
    canonical_edges.reserve(edges.size());
    for (size_t k = 0; k < m*(size_t)n; ++k)
    {
      setword neighbors = cg1[k]; // bit at position l is 1 if there is an edge from k to l
      unsigned int neighbor = 0;
      while (neighbors)
      {
        TAKEBIT(neighbor, neighbors); // get the position of the first 1-bit, and set that bit to 0
        if (k < neighbor) { // avoid duplicates
          canonical_edges.push_back({(uint8_t)k,(uint8_t)neighbor});
        }
      }
    }
    
    // Clean up
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    DYNFREE(map,map_sz);
    DYNFREE(g1,g1_sz);
    DYNFREE(cg1,cg1_sz);

    return true;
  }
}


