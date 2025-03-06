#include "src/disk_quadrangulation.h"
#include "src/cpp_utils.hpp"
#include "src/logging.h"
#include "src/geolog.h"
#include "src/nauty_wrapper.h"
// #include "omp.h"
#include <queue>
#include <assert.h>
#include <bits/stdc++.h>
#include <regex>
#include "gmsh.h"

constexpr bool VERBOSE = false;

/* - Loops */
#define F(_VAR,_NB) for(size_t _VAR = 0; _VAR < (size_t) _NB; ++_VAR)
#define FC(_VAR,_NB,_COND) for(size_t _VAR = 0; _VAR < (size_t) _NB; ++_VAR) if (_COND)
#define RF(...) do {Logging::error(__VA_ARGS__); return false; } while(0)
#define RFC(_COND,...) do { if (_COND) {Logging::error(__VA_ARGS__); return false;} } while(0)

namespace DiskQuadrangulation {
  using std::vector;
  using std::array;
  using std::unordered_map;
  using namespace Logging;

  using Signature = std::array<id,16>;

  Signature quad_mesh_signature(const QuadMesh& M) {
    // TODO: bdr valence + interior valence in two seperate arrays ?
    Signature valence_histogram;
    valence_histogram.fill(0);

    size_t nbv = 0;
    F(i,M.size()) F(j,4) nbv = std::max(nbv,(size_t)M[i][j]+1);
    std::vector<size_t> vert_valence(nbv,0);
    F(i,M.size()) F(j,4) {
      vert_valence[M[i][j]] += 1;
    }

    F(i,vert_valence.size()) {
      id val = vert_valence[i];
      assert(val < 8);
      valence_histogram[val] += 1;
    }

    valence_histogram[0] = M.size(); /* nb quad in first value */

    return valence_histogram;
  }


  inline id2 sorted(id v1, id v2) {
    if (v1 < v2) {
      return {v1,v2};
    } else {
      return {v2,v1};
    }
  }

  inline size_t hash(size_t nb_vert_max, id2 sorted_edge) {
    assert(sorted_edge[0] < nb_vert_max && sorted_edge[1] < nb_vert_max);
    return sorted_edge[0] * nb_vert_max + sorted_edge[1];
  }

  inline size_t hash(size_t nb_vert_max, id4 sorted_quad) {
    assert(sorted_quad[0] < nb_vert_max && sorted_quad[1] < nb_vert_max);
    assert(sorted_quad[2] < nb_vert_max && sorted_quad[3] < nb_vert_max);
    const size_t h = sorted_quad[0] * std::pow(nb_vert_max,3) 
                   + sorted_quad[1] * std::pow(nb_vert_max,2) 
                   + sorted_quad[2] * nb_vert_max 
                   + sorted_quad[3];
    return h;
  }

  inline id4 hash_to_quad(size_t nb_vert_max, size_t hash) {
    const id v4 = hash % nb_vert_max;
    const id v3 = (hash / nb_vert_max) % nb_vert_max;
    const id v2 = (hash / (nb_vert_max*nb_vert_max)) % nb_vert_max;
    const id v1 = (hash / (nb_vert_max*nb_vert_max*nb_vert_max)) % nb_vert_max;
    const id4 quad = {v1,v2,v3,v4};
    return quad;

  }

  vector<size_t> quadMeshHash(size_t nb_vert_max, const QuadMesh& Q) {
    vector<size_t> qhash(Q.size());
    F(i,Q.size()) qhash[i] = hash(nb_vert_max,Q[i]);
    std::sort(qhash.begin(),qhash.end());
    return qhash;
  }

  QuadMesh quad_mesh_hash_to_quad_mesh(size_t nb_vert_max, const vector<size_t>& Q_hash) {
    QuadMesh Q(Q_hash.size());
    F(i,Q.size()) Q[i] = hash_to_quad(nb_vert_max,Q_hash[i]);
    return Q;
  }

  bool is_compatible(size_t nb_vert_max, const std::vector<bool>& Q_edges, const std::vector<bool>& Q_diags, const id4& quad) {
    /* Notes: use hashmap if too slow */
    const size_t diag1 = hash(nb_vert_max,sorted(quad[0],quad[2]));
    const size_t diag2 = hash(nb_vert_max,sorted(quad[1],quad[3]));
    const size_t edge1 = hash(nb_vert_max,sorted(quad[0],quad[1]));
    const size_t edge2 = hash(nb_vert_max,sorted(quad[1],quad[2]));
    const size_t edge3 = hash(nb_vert_max,sorted(quad[2],quad[3]));
    const size_t edge4 = hash(nb_vert_max,sorted(quad[3],quad[0]));
    if (   Q_edges[diag1] || Q_edges[diag2]
        || Q_diags[edge1] || Q_diags[edge2] || Q_diags[edge3] || Q_diags[edge4]
        || Q_diags[diag1] || Q_diags[diag2]
        ) {
      return false;
    }
    return true;
  }

  id4 reorder_quad(id4 quad) {
    id vMin = quad[0];
    id lvMin = 0;
    for (size_t j = 1; j < 4; ++j) {
      if (quad[j] < vMin) {
        vMin = quad[j];
        lvMin = j;
      }
    }
    if (lvMin != 0) {
      id4 quad2 = {
        quad[lvMin],
        quad[(lvMin+1)%4],
        quad[(lvMin+2)%4],
        quad[(lvMin+3)%4]};
      return quad2;
    }
    return quad;
  }

  bool search_shellable_quad_meshes(
      const vector<id>& disk, 
      const vector<size_t>& Q_hash, 
      size_t qMax, size_t nb_vert_max, size_t limit, 
      bool symmetry,
      vector<bool>& already_explored, 
      vector<size_t>& branch_explored, 
      vector<vector<size_t>>& quadMesh_hashes) {

    if (quadMesh_hashes.size() >= limit) return false;
    if ((disk.size()-2)/2 + Q_hash.size() == qMax) return false;
    if (Q_hash.size() == qMax) return false;


    /* Quad edges and diagonals, used to check compatibility */
    id nv = 0;
    F(i,disk.size()) nv = std::max(nv,id(disk[i]+1));
    vector<bool> Q_edges(std::pow(nb_vert_max,2),false);
    vector<bool> Q_diags(std::pow(nb_vert_max,2),false);
    F(f,Q_hash.size()) {
      id4 quad = hash_to_quad(nb_vert_max,Q_hash[f]);
      Q_edges[hash(nb_vert_max,sorted(quad[0],quad[1]))] = true;
      Q_edges[hash(nb_vert_max,sorted(quad[1],quad[2]))] = true;
      Q_edges[hash(nb_vert_max,sorted(quad[2],quad[3]))] = true;
      Q_edges[hash(nb_vert_max,sorted(quad[3],quad[0]))] = true;
      Q_diags[hash(nb_vert_max,sorted(quad[0],quad[2]))] = true;
      Q_diags[hash(nb_vert_max,sorted(quad[1],quad[3]))] = true;
      nv = std::max(nv,id(quad[0]+1));
      nv = std::max(nv,id(quad[1]+1));
      nv = std::max(nv,id(quad[2]+1));
      nv = std::max(nv,id(quad[3]+1));
    }

    if (VERBOSE) {
      vector<int> i_disk(disk.size());
      F(k,i_disk.size()) i_disk[k] = disk[k];
      info("search_shellable_quad_meshes");
      info("  disk = {}", i_disk);
      // info("  Q = {}", Q);
      info("  nv = {}", nv);
      info("  Q_diags = {}", Q_diags);
      info("  Q_edges = {}", Q_edges);
    }

    /* Check if the remaining boundary is a quad */
    if (disk.size() == 4) {
      id4 quad = {disk[0],disk[1],disk[2],disk[3]};
      if (is_compatible(nb_vert_max, Q_edges, Q_diags, quad)) {
        vector<size_t> Qf_hash = Q_hash;
        Qf_hash.push_back(hash(nb_vert_max,reorder_quad(quad)));
        cpp_utils::sort_unique(Qf_hash);
        if (std::find(quadMesh_hashes.begin(),quadMesh_hashes.end(),Qf_hash) == quadMesh_hashes.end()) {
          quadMesh_hashes.push_back(Qf_hash);
        }
      }
    }

    /* Recursive call over compatible flips */
    /*  loop over disk edges */
    F(i,disk.size()) {
      id v1 = disk[i];
      id v2 = disk[(i+1)%disk.size()];

      /*  loop over flip types */
      F(flip_type,3) {
        vector<id> new_disk = disk;
        id4 new_quad;
        if (flip_type == 0) { 
          /* Type 0: create two internal vertices */
          if (nv+2 > nb_vert_max) continue; /* nb_vert_max test */
          id v3 = nv;
          id v4 = nv+1;
          new_quad = {v1,v2,v3,v4};
          if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
          new_disk.insert(new_disk.begin()+i+1,v3);
          new_disk.insert(new_disk.begin()+i+1,v4);
          if (VERBOSE) {
            info("  i={}, flip +2v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
          }
        } else if (flip_type == 1) {
          /* Type 1: create one internal vertex */
          if (nv+1 > nb_vert_max) continue; /* nb_vert_max test */
          id v3 = nv;
          id vPrev = disk[(disk.size()+i-1)%disk.size()];
          new_quad = {v1,v2,v3,vPrev};
          if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
          new_disk[i] = v3;
          if (VERBOSE) {
            info("  i={}, flip +1v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
          }
        } else if (flip_type == 2) {
          /* Type 2: cut the boundary */
          id vPrev = disk[(disk.size()+i-1)%disk.size()];
          id vNext = disk[(i+2)%disk.size()];
          new_quad = {vPrev, v1, v2, vNext};
          if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
          if (i < disk.size()-1) {
            new_disk.erase(new_disk.begin() + i, new_disk.begin() + i + 2);
          } else {
            new_disk.erase(new_disk.begin() + i); /* last vertex of the loop */
            new_disk.erase(new_disk.begin()); /* first vertex of the loop */
          }
          if (VERBOSE) {
            info("  i={}, flip cut, vPrev={}, v1={}, v2={}, vNext={}, new_quad={}", i, vPrev, v1, v2, vNext, new_quad);
          }
        } else {
          continue;
        }
        new_quad = reorder_quad(new_quad);

        if (symmetry && already_explored[hash(nb_vert_max,new_quad)])  continue;

        /* Recursive call */
        vector<size_t> new_Q_hash = Q_hash;
        new_Q_hash.push_back(hash(nb_vert_max,reorder_quad(new_quad)));
        cpp_utils::sort_unique(new_Q_hash);


        vector<size_t> explored_down;
        search_shellable_quad_meshes(new_disk, new_Q_hash, qMax, nb_vert_max, limit,
            symmetry, already_explored, explored_down, quadMesh_hashes);

        if (symmetry) {
          branch_explored.push_back(hash(nb_vert_max,new_quad));

          F(k,explored_down.size()) {
            already_explored[explored_down[k]] = false;
          }

          already_explored[hash(nb_vert_max,new_quad)] = true;
        }
      }
    }

    return true;
  }

  bool find_quad_meshes(
      const std::vector<id>& disk,
      const DiskQuadrangulationOpt& opt,
      std::vector<QuadMesh>& quadMeshes) {

    vector<size_t> Q_empty;
    size_t qmax = opt.nb_quad_max;
    size_t nb_vert_max = opt.nb_vert_max;
    size_t limit = opt.nb_qdrl_max;
    bool symmetry = opt.symmetry;

    if (nb_vert_max >= ID_MAX) {
      error("cannot deal with more than {} vertices", int(ID_MAX));
      return false;
    }

    const size_t hash_max = std::pow(nb_vert_max,4);
    if (symmetry) {
      info("hash_max_id4 = {}", hash_max);
    }
    info("hash_max_id2 = {}", nb_vert_max*nb_vert_max);
    vector<bool> already_explored(hash_max,false);
    vector<size_t> branch_explored;
    vector<vector<size_t>> quadMesh_hashes;
    search_shellable_quad_meshes(disk, Q_empty, qmax, nb_vert_max, limit, symmetry, already_explored, branch_explored, quadMesh_hashes);

    quadMeshes.resize(quadMesh_hashes.size());
    F(i,quadMesh_hashes.size()) {
      quadMeshes[i] = quad_mesh_hash_to_quad_mesh(nb_vert_max,quadMesh_hashes[i]);
    }

    return true;
  }

  using QuadMeshHash = std::vector<size_t> ;

  bool generate_quad_meshes(
      const DiskQuadrangulationOpt& opt,
      std::vector<QuadMesh>& quadMeshes) {
    info("generate quad meshes ... nb quad max = {}", opt.nb_quad_max);

    const size_t nb_vert_max = opt.nb_vert_max;
    info("nb_vert_max = {}, ^4 = {} (used for hashing)", nb_vert_max, std::pow(nb_vert_max,4));

    /* TODO: use fixed size vector to store the quads in a mesh ? qmax known */

    vector<QuadMeshHash> meshes;
    vector<Signature> signatures;
    std::queue<std::pair<QuadMeshHash,vector<id>>> queue;


    id4 quad0 = {0,1,2,3};
    QuadMeshHash Q0;
    Q0.push_back(hash(nb_vert_max,reorder_quad(quad0)));
    Signature sig0 = quad_mesh_signature(quad_mesh_hash_to_quad_mesh(nb_vert_max,Q0));
    signatures.push_back(sig0);
    vector<id> bdrLoop0 = {0,1,2,3};
    queue.push({Q0,bdrLoop0});
    while (!queue.empty()) {
      QuadMeshHash Q = queue.front().first;
      vector<id> bdrLoop = queue.front().second;
      queue.pop();
      meshes.push_back(Q);

      id4 quad = hash_to_quad(nb_vert_max,Q[0]);

      if (Q.size() == opt.nb_quad_max) continue;
      id nv = 0;
      vector<bool> Q_edges(std::pow(nb_vert_max,2),false);
      vector<bool> Q_diags(std::pow(nb_vert_max,2),false);
      F(f,Q.size()) {
        id4 quad = hash_to_quad(nb_vert_max,Q[f]);
        Q_edges[hash(nb_vert_max,sorted(quad[0],quad[1]))] = true;
        Q_edges[hash(nb_vert_max,sorted(quad[1],quad[2]))] = true;
        Q_edges[hash(nb_vert_max,sorted(quad[2],quad[3]))] = true;
        Q_edges[hash(nb_vert_max,sorted(quad[3],quad[0]))] = true;
        Q_diags[hash(nb_vert_max,sorted(quad[0],quad[2]))] = true;
        Q_diags[hash(nb_vert_max,sorted(quad[1],quad[3]))] = true;
        nv = std::max(nv,id(quad[0]+1));
        nv = std::max(nv,id(quad[1]+1));
        nv = std::max(nv,id(quad[2]+1));
        nv = std::max(nv,id(quad[3]+1));
      }

      /* Loop over flips */
      F(i,bdrLoop.size()) {
        id v1 = bdrLoop[i];
        id v2 = bdrLoop[(i+1)%bdrLoop.size()];

        /*  loop over flip types */
        F(flip_type,3) {
          vector<id> new_bdrLoop = bdrLoop;
          id4 new_quad;
          if (flip_type == 0) { 
            /* Type 0: create two external vertices */
            id v3 = nv;
            id v4 = nv+1;
            new_quad = {v1,v2,v3,v4};
            if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
            new_bdrLoop.insert(new_bdrLoop.begin()+i+1,v3);
            new_bdrLoop.insert(new_bdrLoop.begin()+i+1,v4);
            if (VERBOSE) {
              info("  i={}, flip +2v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
            }
          } else if (flip_type == 1) {
            /* Type 1: create one internal vertex */
            id v3 = nv;
            id vPrev = bdrLoop[(bdrLoop.size()+i-1)%bdrLoop.size()];
            new_quad = {v1,v2,v3,vPrev};
            if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
            new_bdrLoop[i] = v3;
            if (VERBOSE) {
              info("  i={}, flip +1v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
            }
          } else if (flip_type == 2) {
            /* Type 2: cut the boundary */
            id vPrev = bdrLoop[(bdrLoop.size()+i-1)%bdrLoop.size()];
            id vNext = bdrLoop[(i+2)%bdrLoop.size()];
            new_quad = {vPrev, v1, v2, vNext};
            if (!is_compatible(nb_vert_max,Q_edges,Q_diags,new_quad)) continue; /* compat. test */
            if (i < bdrLoop.size()-1) {
              new_bdrLoop.erase(new_bdrLoop.begin() + i, new_bdrLoop.begin() + i + 2);
            } else {
              new_bdrLoop.erase(new_bdrLoop.begin() + i); /* last vertex of the loop */
              new_bdrLoop.erase(new_bdrLoop.begin()); /* first vertex of the loop */
            }
            if (VERBOSE) {
              info("  i={}, flip cut, vPrev={}, v1={}, v2={}, vNext={}, new_quad={}", i, vPrev, v1, v2, vNext, new_quad);
            }
          } else {
            continue;
          }
          new_quad = reorder_quad(new_quad);
          QuadMeshHash new_Q = Q;
          new_Q.push_back(hash(nb_vert_max,reorder_quad(new_quad)));
          cpp_utils::sort_unique(new_Q);

          Signature candidate = quad_mesh_signature(quad_mesh_hash_to_quad_mesh(nb_vert_max,new_Q));
          if (cpp_utils::inVector(candidate, signatures)) continue;

          queue.push({new_Q,new_bdrLoop});
          signatures.push_back(candidate);
        }
      }
    }

    quadMeshes.resize(meshes.size());
    F(i,meshes.size()) {
      quadMeshes[i] = quad_mesh_hash_to_quad_mesh(nb_vert_max,meshes[i]);
    }

    return true;
  }

  struct id2Hash {
    size_t operator()(id2 p) const noexcept {
      return size_t(p[0]) << 32 | p[1];
    }
  };

  struct vidHash {
    size_t operator()(const std::vector<id>& p) const noexcept {
      uint32_t hash = 0;
      for (size_t i = 0; i < p.size(); ++i) {
        hash += p[i];
        hash += hash << 10;
        hash ^= hash >> 6;
      }
      hash += hash << 3;
      hash ^= hash >> 11;
      hash += hash << 15;
      return hash;
    }
  };

  bool get_smallest_rotation(const vector<id>& bdrVal, const vector<id>& bdrLoop, 
      vector<id>& smallestBdrVal, vector<id>& smallestBdrLoop) {
    smallestBdrVal = bdrVal;
    smallestBdrLoop = bdrLoop;
    vector<id> rot = bdrVal;
    vector<id> ids = bdrLoop;
    F(i,bdrVal.size()) {
      std::rotate(rot.begin(),rot.begin()+1,rot.end());
      std::rotate(ids.begin(),ids.begin()+1,ids.end());
      if (rot < smallestBdrVal) {
        smallestBdrVal = rot;
        smallestBdrLoop = ids;
      }
    }
    std::reverse(rot.begin(),rot.end());
    std::reverse(ids.begin(),ids.end());
    F(i,bdrVal.size()) {
      std::rotate(rot.begin(),rot.begin()+1,rot.end());
      std::rotate(ids.begin(),ids.begin()+1,ids.end());
      if (rot < smallestBdrVal) {
        smallestBdrVal = rot;
        smallestBdrLoop = ids;
      }
    }
    return true;
  }

  struct QMesh {
    id n;
    vector<id4> quads;
    vector<id> loop;
  };

  bool makeValenceCanonical(QMesh& M) {
    std::vector<id> qval(M.n,0);
    F(i,M.quads.size()) F(j,4) qval[M.quads[i][j]] += 1;

    vector<id> bdr_val_loop(M.loop.size());
    F(i,M.loop.size()) bdr_val_loop[i] = qval[M.loop[i]];

    /* Smallest 'boundary val' loop */
    vector<id> s_bdrVal(bdr_val_loop.size());
    vector<id> s_overt(bdr_val_loop.size());
    get_smallest_rotation(bdr_val_loop, M.loop, s_bdrVal, s_overt);

    /* Re-numbering map */
    vector<id> old2new(M.n,NO_ID);
    id cv = 0;
    F(i,M.loop.size()) {
      old2new[M.loop[i]] = cv;
      cv += 1;
    }
    FC(i,M.n,old2new[i] == NO_ID) {
      old2new[i] = cv;
      cv += 1;
    }

    /* Apply */
    F(i,M.quads.size()) F(j,4) M.quads[i][j] = old2new[M.quads[i][j]];
    F(i,M.loop.size()) M.loop[i] = old2new[M.loop[i]];

    return true;
  }

  std::ostream& operator<<(std::ostream& os, const QMesh& M) { 
    os << "(";
    os << "n:" << (uint32_t) M.n;
    os << ",quads:" << M.quads;
    os << ",bdrLoop:" << M.loop;
    os << ")";
    return os;
  }


  void qmesh_valences(const QMesh& M, vector<id>& int_val_nb, vector<id>& bdr_val_nb) {
    int_val_nb.clear();
    bdr_val_nb.clear();
    int_val_nb.resize(5+M.quads.size(),0);
    bdr_val_nb.resize(3+M.quads.size(),0);
    std::vector<bool> onBdr(M.n,false);
    F(i,M.loop.size()) onBdr[M.loop[i]] = true;
    std::vector<id> qval(M.n,0);
    F(i,M.quads.size()) F(j,4) qval[M.quads[i][j]] += 1;
    id nb = 0;
    F(v,M.n) {
      if (onBdr[v]) {
        bdr_val_nb[qval[v]] += 1;
        nb += 1;
      } else {
        int_val_nb[qval[v]] += 1;
      }
    }
    return;
  }

  void qmesh_max_valences(const QMesh& M, id& int_val_max, id& bdr_val_max) {
    int_val_max = 0;
    bdr_val_max = 0;
    std::vector<id> int_val_nb;
    std::vector<id> bdr_val_nb;
    qmesh_valences(M, int_val_nb, bdr_val_nb);
    FC(k,int_val_nb.size(), int_val_nb[k] > 0) {
      int_val_max = std::max(int_val_max,id(k));
    }
    FC(k,bdr_val_nb.size(), bdr_val_nb[k] > 0) {
      bdr_val_max = std::max(bdr_val_max,id(k));
    }
  }


  /* format: [nb_int_vert, valence1, #, valence2, #, ..., nb_bdr_vert, bvalence1, #, bvalence2, #, ...] */
  void qmesh_signature(const QMesh& M, std::vector<id>& sig) {
    std::vector<id> int_val_nb;
    std::vector<id> bdr_val_nb;
    qmesh_valences(M, int_val_nb, bdr_val_nb);
    sig.clear();
    sig.reserve(12);
    sig.push_back(M.n);
    FC(i,int_val_nb.size(),int_val_nb[i] > 0) {
      sig.push_back(i);
      sig.push_back(int_val_nb[i]);
    }
    size_t nb = M.loop.size();
    sig.push_back(nb);
    FC(i,bdr_val_nb.size(),bdr_val_nb[i] > 0) {
      sig.push_back(i);
      sig.push_back(bdr_val_nb[i]);
    }
  }

  using vec3 = std::array<double,3>;
  inline vec3 operator-(const vec3& a, const vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
  inline vec3 operator+(const vec3& a, const vec3& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
  inline vec3 operator*(const double& a, const vec3& b) { return {a*b[0], a*b[1], a*b[2]}; }
  inline vec3 operator*(const vec3& a, const double& b) { return {a[0]*b, a[1]*b, a[2]*b}; }

  bool laplacian_smoothing(const QMesh& M, vector<vec3>& points, size_t iter, double lambda) {
    vector<bool> locked(M.n);
    for (id v: M.loop) locked[v] = true;
    vector<vector<id>> v2v(M.n);
    F(f,M.quads.size()) F(le,4) {
      v2v[M.quads[f][le]].push_back(M.quads[f][(le+1)%4]);
      v2v[M.quads[f][(le+1)%4]].push_back(M.quads[f][le]);
    }
    F(v,M.n) cpp_utils::sort_unique(v2v[v]);

    F(k,iter) {
      FC(v,v2v.size(),v2v[v].size() > 0 && !locked[v]) {
        vec3 avg = {0.,0.,0.};
        F(lv,v2v[v].size()) {
          id v2 = v2v[v][lv];
          vec3 p2 = points[v2];
          avg = avg + p2;
        }
        vec3 target = 1./v2v[v].size() * avg;
        points[v] = lambda * target + (1.-lambda) * points[v];
      }
    }
    return true;
  }

  bool visualize_quadrangulation(vec3 p, double r, const QMesh& M, const std::string& view, bool showCorners = true) {
    vector<vec3> points(M.n);
    vector<double> colors(M.n,0);

    if (showCorners) {
      std::vector<id> qval(M.n,0);
      F(i,M.quads.size()) F(j,4) qval[M.quads[i][j]] += 1;
      F(i,M.loop.size()) {
        id v = M.loop[i];
        if (qval[v] == 1) colors[v] = 1.;
      }
    }

    F(i,M.loop.size()) {
      id v = M.loop[i];
      double agl = 2. * M_PI * double(i) / double(M.loop.size());
      vec3 p2 = {p[0] + r * cos(agl), p[1] + r * sin(agl), 0};
      if (showCorners && colors[v] == 0) {
        double factor = std::sqrt(2)/2;
        p2 = {p[0] + factor * r * cos(agl), p[1] + factor * r * sin(agl), 0};
      }
      points[v] = p2;
    }
    laplacian_smoothing(M, points, 100, 1.);
    F(f,M.quads.size()) {
      vector<vec3> pts = {points[M.quads[f][0]], points[M.quads[f][1]],
        points[M.quads[f][2]], points[M.quads[f][3]] };
      vector<double> values = {colors[M.quads[f][0]], colors[M.quads[f][1]],
        colors[M.quads[f][2]], colors[M.quads[f][3]] };
      GeoLog::add(pts,values,view);
    }
    return true;
  }

  bool visualize_quadrangulations(const std::vector<QMesh>& meshes, bool filterConcave = false) {
    std::map<id2 ,vector<uint32_t>> nb_ni_to_mids;
    F(i,meshes.size()) {
      const QMesh& M = meshes[i];
      if (filterConcave) {
        id int_val_max = 0;
        id bdr_val_max = 0;
        qmesh_max_valences(M, int_val_max, bdr_val_max);
        if (bdr_val_max > 2) continue;
      }
      id nb = M.loop.size();
      id ni = M.n - M.loop.size();
      id2 vp = {nb,ni};
      nb_ni_to_mids[vp].push_back(i);
    }

    double x = 0.;
    double y = 0.;
    for (const auto& kv: nb_ni_to_mids) {
      std::string name = fmt::format("{}bdr, {}int",kv.first[0],kv.first[1]);
      x = 0.;
      y -= 1.;
      trace("viz | view {} at y={} with {} meshes",name,y,kv.second.size());
      F(i,kv.second.size()) {
        uint32_t mid = kv.second[i];
        const QMesh& M = meshes[mid];
        id nb = M.loop.size();
        id ni = M.n - M.loop.size();
        if (nb != kv.first[0] || ni != kv.first[1]) {
          error("bad ! kv.first={} but nb,ni={},{}",kv.first,nb,ni);
          return false;
        }
        x = double(i);
        visualize_quadrangulation({x,y,0},0.45,M,name);
      }
      GeoLog::flush();
    }
    GeoLog::flush();
    vector<int> tags;
    gmsh::view::getTags(tags);
    for (id tag: tags) {
      gmsh::option::setNumber("View["+std::to_string(tag)+"].Visible", 1);
      // gmsh::option::setNumber("View["+std::to_string(tag)+"].Explode", 1.);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].ShowElement", 1);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].DrawSkinOnly", 1);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].LineWidth", 1.);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].LineType", 1);
    }

    return true;
  }

  bool visualize_concave_quadrangulations(const std::vector<QMesh>& meshes) {
    std::map<id ,vector<uint32_t>> nc_to_mids;
    F(i,meshes.size()) {
      const QMesh& M = meshes[i];
      id int_val_max = 0;
      id bdr_val_max = 0;
      qmesh_max_valences(M, int_val_max, bdr_val_max);
      if (bdr_val_max > 2) continue; /* not concave */

      std::vector<id> int_val_nb;
      std::vector<id> bdr_val_nb;
      qmesh_valences(M, int_val_nb, bdr_val_nb);
      id nc = bdr_val_nb[1];
      nc_to_mids[nc].push_back(i);
    }

    double x = 0.;
    double y = 0.;
    for (const auto& kv: nc_to_mids) {
      std::string name = fmt::format("{} corners",kv.first);
      x = 0.;
      y -= 1.;
      trace("viz | view {} at y={} with {} meshes",name,y,kv.second.size());
      F(i,kv.second.size()) {
        uint32_t mid = kv.second[i];
        const QMesh& M = meshes[mid];
        id nb = M.loop.size();
        id ni = M.n - M.loop.size();
        x = double(i);
        visualize_quadrangulation({x,y,0},0.45,M,name);
        GeoLog::add({x,y,0}, double(mid+1), name + "_numbers");
      }
      GeoLog::flush();
    }
    GeoLog::flush();
    vector<int> tags;
    gmsh::view::getTags(tags);
    for (id tag: tags) {
      gmsh::option::setNumber("View["+std::to_string(tag)+"].Visible", 1);
      // gmsh::option::setNumber("View["+std::to_string(tag)+"].Explode", 1.);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].ShowElement", 1);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].DrawSkinOnly", 1);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].LineWidth", 1.);
      gmsh::option::setNumber("View["+std::to_string(tag)+"].LineType", 1);
    }

    return true;
  }

  std::string ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
  }

  // bool export_remeshing_patterns(const std::vector<QMesh>& meshes, const std::string& path) {
  //   std::map<vector<id> ,vector<QMesh> > bdrValLoop_to_mids;

  //   F(i,meshes.size()) {
  //     QMesh M = meshes[i]; /* copy */

  //     /* Get valences on bdr loop, rotated with minimal values first */
  //     makeValenceCanonical(M);
  //     std::vector<id> qval(M.n,0);
  //     F(i,M.quads.size()) F(j,4) qval[M.quads[i][j]] += 1;
  //     vector<id> bdr_val_loop(M.loop.size());
  //     F(i,M.loop.size()) bdr_val_loop[i] = qval[M.loop[i]];

  //     bdrValLoop_to_mids[bdr_val_loop].push_back(M);
  //   }

  //   DBG(bdrValLoop_to_mids.size());
  //   for (const auto& kv: bdrValLoop_to_mids) {
  //     const vector<id>& bdr_val_loop = kv.first;
  //     if (kv.second.size() > 1) {
  //       DBG(bdr_val_loop.size(),bdr_val_loop, kv.second.size());
  //     }
  //   }

  //   return true;
  // }

  bool export_quadrangulations(const std::vector<QMesh>& meshes, const std::string& path, bool cpp = false, bool manual = false) {
    info("export {} quadrangulations to {}", meshes.size(),path);
    std::map<id2 ,vector<uint32_t>> nb_ni_to_mids;
    F(i,meshes.size()) {
      const QMesh& M = meshes[i];
      id nb = M.loop.size();
      id ni = M.n - M.loop.size();
      id2 vp = {nb,ni};
      nb_ni_to_mids[vp].push_back(i);
    }
    std::ofstream out;
    out.open (path);
    if (cpp) {
      for (const auto& kv: nb_ni_to_mids) {
        id B = kv.first[0];
        id I = kv.first[1];
        std::string content = "const std::vector< std::vector< std::array<uint32_t,4> > > disk_quadrangulations_B"+std::to_string(B) + "_I" + std::to_string(I) + " = {\n";
        F(i,kv.second.size()) {
          uint32_t mid = kv.second[i];
          const QMesh& M = meshes[mid];
          id nb = M.loop.size();
          id ni = M.n - M.loop.size();
          if (nb != kv.first[0] || ni != kv.first[1]) {
            error("bad ! kv.first={} but nb,ni={},{}",kv.first,nb,ni);
            return false;
          }
          std::string str = fmt::format("{}",M.quads);
          str = ReplaceAll(str, "[", "{");
          str = ReplaceAll(str, "]", "}");
          content += "    " + str + ",\n";
        }
        content += "};\n";
        out << content;
      }
    } else{
      bool big_string_literal = false;
      if (big_string_literal) {
        out << "// #v_boundary #v_interior #quads q1v1 q1v2 q1v3 q1v4 q2v1 q2v2 q2v3 q2v4 ... qNv1 qNv2 qNv3 qNv4\n";
        for (const auto& kv: nb_ni_to_mids) {
          id B = kv.first[0];
          id I = kv.first[1];
          F(i,kv.second.size()) {
            uint32_t mid = kv.second[i];
            const QMesh& M = meshes[mid];
            id nb = M.loop.size();
            id ni = M.n - M.loop.size();
            if (nb != kv.first[0] || ni != kv.first[1]) {
              error("bad ! kv.first={} but nb,ni={},{}",kv.first,nb,ni);
              return false;
            }
            out << (uint32_t)B << " " << (uint32_t)I << " " << (uint32_t)M.quads.size();
            F(f,M.quads.size()) F(lv,4) {
              out << " " <<  (uint32_t)M.quads[f][lv];
            }
            out << "\n";
          }
        }
      } else {
        const size_t SPLIT_FOR_MSVC = 500;
        out << "#pragma once\n";
        out << "#include <string>\n\n";
        out << "// #v_boundary #v_interior #quads q1v1 q1v2 q1v3 q1v4 q2v1 q2v2 q2v3 q2v4 ... qNv1 qNv2 qNv3 qNv4\n";
        size_t count = 0;
        size_t nblock = 0;
        for (const auto& kv: nb_ni_to_mids) {
          id B = kv.first[0];
          id I = kv.first[1];
          F(i,kv.second.size()) {
            uint32_t mid = kv.second[i];
            const QMesh& M = meshes[mid];
            id nb = M.loop.size();
            id ni = M.n - M.loop.size();
            if (nb != kv.first[0] || ni != kv.first[1]) {
              error("bad ! kv.first={} but nb,ni={},{}",kv.first,nb,ni);
              return false;
            }

            if (manual) {
              std::vector<id> int_val_nb;
              std::vector<id> bdr_val_nb;
              qmesh_valences(M, int_val_nb, bdr_val_nb);

              /* If only irregular vertices inside, avoid */
              if (nb >= 10 && ni >= 3 && int_val_nb[4] == 0) continue;
              if (nb >= 18 && ni >= 2) continue;
              if (nb >= 20 && ni >= 1) continue;
            }




            if (count % SPLIT_FOR_MSVC == 0) {
              if (nblock > 0) { /* close current block */
                out << ";\n\n";
              }
              /* Start a new const char* */
              out << "const char* data_dqrgl_block" + std::to_string(nblock) + " = \n";
              nblock += 1;
            }

            out << "  \"";
            out << (uint32_t)B << " " << (uint32_t)I << " " << (uint32_t)M.quads.size();
            F(f,M.quads.size()) F(lv,4) {
              out << " " <<  (uint32_t)M.quads[f][lv];
            }
            out << "\\n\"\n";
            count += 1;
          }
        }
        if (nblock > 0) { /* close current block */
          out << ";\n\n";
        }
        std::string concat_func = "void diskQuadrangulationConcat(std::string& data) {\n";
        concat_func += "  data = data ";
        F(i,nblock) {
          concat_func += "+ std::string(data_dqrgl_block" + std::to_string(i) + ") ";
        }
        concat_func += ";\n};\n";
        out << concat_func << "\n";
      }
    }
    out.close();

    return true;
  }

  std::string string_valence_repartition(const QMesh& M) {
    std::vector<id> int_val_nb;
    std::vector<id> bdr_val_nb;
    qmesh_valences(M, int_val_nb, bdr_val_nb);
    std::string msg = fmt::format("{}v {}q | boundary: {} vert. (#2={}",
        M.n,M.quads.size(),M.loop.size(),bdr_val_nb[2]);
    FC(i,bdr_val_nb.size(),i!=2&&bdr_val_nb[i]>0) {
      msg += fmt::format(", #{}={}", i, bdr_val_nb[i]);
    }
    msg += fmt::format(") | interior: {} vert.",M.n-M.loop.size());
    if (M.n > M.loop.size()) {
      msg += fmt::format(" (#4={}",int_val_nb[4]);
      FC(i,int_val_nb.size(),i!=4&&int_val_nb[i]>0) {
        msg += fmt::format(", #{}={}", i, int_val_nb[i]);
      }
      msg += ")";
    }
    return msg;
  }

  inline size_t edge_no(id v1, id v2) {
    if (v1 < v2) {
      return size_t(v1) * size_t(ID_MAX) + size_t(v2);
    }
    return size_t(v2) * size_t(ID_MAX) + size_t(v1);
  }

  constexpr size_t ID_MAX_2 = ID_MAX * ID_MAX;

  bool is_compatible(const std::bitset<ID_MAX_2>& Q_edges, const std::bitset<ID_MAX_2>& Q_diags, 
      const id4& quad) {
    const size_t diag1 = edge_no(quad[0],quad[2]);
    const size_t diag2 = edge_no(quad[1],quad[3]);
    const size_t edge1 = edge_no(quad[0],quad[1]);
    const size_t edge2 = edge_no(quad[1],quad[2]);
    const size_t edge3 = edge_no(quad[2],quad[3]);
    const size_t edge4 = edge_no(quad[3],quad[0]);
    if (   Q_edges.test(diag1) || Q_edges.test(diag2)
        || Q_diags.test(edge1) || Q_diags.test(edge2) || Q_diags.test(edge3) || Q_diags.test(edge4)
        || Q_diags.test(diag1) || Q_diags.test(diag2)
        ) {
      return false;
    }
    return true;
  }

  struct QuadMeshStorage {
    std::vector<QMesh> meshes;
    std::vector<std::unordered_set<std::vector<id>,vidHash> > n_signatures;
    std::vector<std::unordered_set<std::vector<id>,vidHash> > n_canonical_edges;

    bool equivalent_already_inside(const QMesh& M, bool nauty_check = true) {
      size_t n = M.n;
      vector<id> sig;
      vector<id> ces;
      if (n < n_signatures.size()) {
        /* Check if existing mesh with same valence signature */
        qmesh_signature(M, sig);
        auto its = n_signatures[n].find(sig);
        if (its != n_signatures[n].end()) {
          if (nauty_check) {
            /* Check if existing mesh with same canonical edges */
            vector<id2> edges(4*M.quads.size());
            F(f,M.quads.size()) F(le,4) {
              edges[4*f+le] = sorted(M.quads[f][le],M.quads[f][(le+1)%4]);
            }
            cpp_utils::sort_unique(edges);
            vector<id2> canonical_edges;
            bool ok = Nauty::graph_canonical_labelling(M.n,edges,canonical_edges);
            RFC(!ok, "failed to get canonical labelling");
            ces.resize(2*canonical_edges.size());
            F(i,canonical_edges.size()) F(j,2) ces[2*i+j] = canonical_edges[i][j];
            auto itc = n_canonical_edges[n].find(ces);
            if (itc != n_canonical_edges[n].end()) {
              /* Equivalent mesh exists ! */
              return true;
            }
          } else {
            return true;
          }
        }
      }
      return false;
    }

    bool add_quad_mesh_if_no_equivalent(const QMesh& M, bool nauty_check = true) {
      size_t n = M.n;
      vector<id> sig;
      vector<id> ces;
      if (n < n_signatures.size()) {
        /* Check if existing mesh with same valence signature */
        qmesh_signature(M, sig);
        auto its = n_signatures[n].find(sig);
        if (its != n_signatures[n].end()) {
          if (nauty_check) {
            /* Check if existing mesh with same canonical edges */
            vector<id2> edges(4*M.quads.size());
            F(f,M.quads.size()) F(le,4) {
              edges[4*f+le] = sorted(M.quads[f][le],M.quads[f][(le+1)%4]);
            }
            cpp_utils::sort_unique(edges);
            vector<id2> canonical_edges;
            bool ok = Nauty::graph_canonical_labelling(M.n,edges,canonical_edges);
            RFC(!ok, "failed to get canonical labelling");
            ces.resize(2*canonical_edges.size());
            F(i,canonical_edges.size()) F(j,2) ces[2*i+j] = canonical_edges[i][j];
            auto itc = n_canonical_edges[n].find(ces);
            if (itc != n_canonical_edges[n].end()) {
              /* Equivalent mesh exists ! */
              trace("do not add bc. mesh with same canonical edges exists");
              return false;
            }
          } else {
            trace("do not add bc. mesh with same valence signature exists and no nauty check");
            return false;
          }
        }
      }

      /* Add the mesh to storage */
      if (sig.size() == 0) {
        qmesh_signature(M, sig);
      }
      if (n >= n_signatures.size()) {
        n_signatures.resize(n+1);
      }
      if (nauty_check) {
        if (ces.size() == 0) {
          vector<id2> edges(4*M.quads.size());
          F(f,M.quads.size()) F(le,4) {
            edges[4*f+le] = sorted(M.quads[f][le],M.quads[f][(le+1)%4]);
          }
          cpp_utils::sort_unique(edges);
          vector<id2> canonical_edges;
          bool ok = Nauty::graph_canonical_labelling(M.n,edges,canonical_edges);
          RFC(!ok, "failed to get canonical labelling");
          ces.resize(2*canonical_edges.size());
          F(i,canonical_edges.size()) F(j,2) ces[2*i+j] = canonical_edges[i][j];
        }
        if (n >= n_canonical_edges.size()) {
          n_canonical_edges.resize(n+1);
        }
        n_canonical_edges[n].insert(ces);
      }
      n_signatures[n].insert(sig);
      meshes.push_back(M);

      return true;
    }

  };

  void process_mesh_stats(const QMesh& M, const DiskQuadrangulationOpt& opt, 
      bool& add_to_meshes, bool& add_to_queue) {
    add_to_meshes = true;
    add_to_queue = true;

    /* Check global stats (number of vertices) */
    id nv_int = M.n - M.loop.size();
    if (M.n > opt.nb_vert_max || nv_int > opt.nb_int_vert_max || M.quads.size() >= opt.nb_quad_max) {
      add_to_meshes = false;
      add_to_queue = false;
    } else if (M.loop.size() > opt.nb_bdr_vert_max) {
      add_to_meshes = false;
      add_to_queue = true; /* flip cut may reduce number of bdr vertices */
    }
    if (!add_to_meshes && !add_to_queue) return;

    /* Check valence restrictions */
    if (opt.val_bdr_max == ID_MAX && opt.val_int_max == ID_MAX) return;
    id int_val_max = 0;
    id bdr_val_max = 0;
    qmesh_max_valences(M, int_val_max, bdr_val_max);
    if (opt.val_int_max != ID_MAX && int_val_max > opt.val_int_max) {
      /* interior valence will not decrease */
      add_to_meshes = false;
      add_to_queue = false;
      return;
    }
    if (opt.val_bdr_max != ID_MAX && bdr_val_max > opt.val_bdr_max) {
      add_to_meshes = false;
      /* boundary valence may decrease (max bdr val becomes interior vertex) */
      add_to_queue = true;
    }

    if (opt.nb_quad_max > 2 && M.quads.size() == opt.nb_quad_max - 1) {
      add_to_queue = false; 
    }

    return;
  }


  bool enumerate_disk_quadrangulations(const DiskQuadrangulationOpt& opt) {
    info("enumerate disk quadrangulations ...");
    if (opt.no_nauty_check) {
      warn("Nauty check (i.e. graph equivalence) disabled, only using valence signature");
    }

    // TODO issue with memory consumption, something is getting too big and not cleared / shrinked ?

    std::vector<QMesh> meshes;
    std::vector<bool> meshes_add_store;
    std::vector<bool> meshes_add_queue;
    std::queue<QMesh> queue;

    /* Storage for unique quad meshes */
    QuadMeshStorage store;
    QuadMeshStorage explored; /* to keep track of explored pathes, already in queue */

    /* Init the queue with a quad */
    QMesh M0;
    M0.n = 4;
    M0.quads = {reorder_quad({0,1,2,3})};
    M0.loop = {0,1,2,3};
    queue.push(M0);
    trace("init queue with M0={}",M0);
    store.add_quad_mesh_if_no_equivalent(M0,!opt.no_nauty_check);

    std::bitset<ID_MAX_2> Q_edges;
    std::bitset<ID_MAX_2> Q_diags;

    vector<bool> combi_reached(ID_MAX_2,false);
    while (!queue.empty()) {
      QMesh M = queue.front();
      queue.pop();

      trace("---");
      trace("queue size: {}, current quadrangulation: N={},B={},Q={}", queue.size(), M.n, M.loop.size(),M.quads.size());

      if (M.quads.size() >= opt.nb_quad_max) {
        trace("reached maximum number of quads: {}", opt.nb_quad_max);
        continue;
      }

      if (opt.verbose) {
        size_t combi = ID_MAX*size_t(M.loop.size())+size_t(M.n-M.loop.size());
        if (!combi_reached[combi]) {
          info("- reaching #boundary={}, #interior={} ({} meshes in queue, {} in explored, {} in store)...",M.loop.size(),M.n - M.loop.size(), 
              queue.size(), explored.meshes.size(), store.meshes.size());
          combi_reached[combi] = true;
        }
      }


      /* Flag quad edges and diagonals of current mesh */
      vector<id> qVal(M.n,0);
      Q_edges.reset();
      Q_diags.reset();
      F(f,M.quads.size()) {
        const id4& quad = M.quads[f];
        Q_edges.set(edge_no(quad[0],quad[1]));
        Q_edges.set(edge_no(quad[1],quad[2]));
        Q_edges.set(edge_no(quad[2],quad[3]));
        Q_edges.set(edge_no(quad[3],quad[0]));
        Q_diags.set(edge_no(quad[0],quad[2]));
        Q_diags.set(edge_no(quad[1],quad[3]));
        F(lv,4) qVal[quad[lv]] += 1;
      }

      /* Loop over flips */
      meshes.clear(); /* store all meshes built from the flips */
      meshes_add_queue.clear();
      meshes_add_store.clear();

      trace("loop over flips, bdrLoop = {} ...", M.loop);
      id nv = M.n;
      F(i,M.loop.size()) {
        id v1 = M.loop[i];
        id v2 = M.loop[(i+1)%M.loop.size()];

        /*  loop over flip types */
        F(flip_type,3) {
          vector<id> new_bdrLoop = M.loop;
          id4 new_quad;
          id n2 = nv;

          if (flip_type == 0) { 
            /* Type 0: create two external vertices, no new internal vertex */
            if (nv+2 > opt.nb_vert_max) {
              trace("- i={}, flip +2v canceled, would be over maximum number of vertices: {}", i, opt.nb_vert_max);
              continue;
            }
            id v3 = nv;
            id v4 = nv+1;
            new_quad = {v2,v1,v4,v3};
            n2 = M.n + 2;
            if (!is_compatible(Q_edges,Q_diags,new_quad)) {
              trace("- i={}, flip +2v canceled, incompatible", i);
              continue; /* compat. test */
            }
            new_bdrLoop.insert(new_bdrLoop.begin()+i+1,v3);
            new_bdrLoop.insert(new_bdrLoop.begin()+i+1,v4);
            trace("- i={}, flip +2v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
          } else if (flip_type == 1) {
            /* Type 1: create one external vertex, one bdr vertex become interior vertex */
            if (nv+1 > opt.nb_vert_max) {
              trace("- i={}, flip +1v canceled, would be over maximum number of vertices: {}", i, opt.nb_vert_max);
              continue;
            }
            id v3 = nv;
            id vPrev = M.loop[(M.loop.size()+i-1)%M.loop.size()];
            // new_quad = {v1,v2,v3,vPrev};
            new_quad = {v2,v1,vPrev,v3};
            n2 = M.n + 1;
            /* v1 becomes interior vertex, check nb and valence */
            if (opt.nb_int_vert_max != ID_MAX && M.n - M.loop.size() + 1 > opt.nb_int_vert_max) {
              /* the flip would increase the interior number of vertices over the max, cancel */
              trace("- i={}, flip +1v canceled, interior nb of vertices would be too high", i);
              continue;
            }
            if (opt.val_int_max != ID_MAX && qVal[v1] >= opt.val_int_max) {
              /* the flip would increase the interior valence over the max, cancel */
              trace("- i={}, flip +1v canceled, interior valence would be too high", i);
              continue;
            }
            if (!is_compatible(Q_edges,Q_diags,new_quad)) {
              trace("- i={}, flip +1v canceled, incompatible", i);
              continue; /* compat. test */
            }
            new_bdrLoop[i] = v3;
            trace("- i={}, flip +1v, v1={}, v2={}, new_quad={}", i, v1, v2, new_quad);
          } else if (flip_type == 2) {
            /* Type 2: cut the boundary, v1, v2 become interior vertices */
            id vPrev = M.loop[(M.loop.size()+i-1)%M.loop.size()];
            id vNext = M.loop[(i+2)%M.loop.size()];
            new_quad = {vNext, v2, v1, vPrev};
            /* v1,v2 become interior vertex, check nb and valence */
            if (opt.nb_int_vert_max != ID_MAX && M.n - M.loop.size() + 2 > opt.nb_int_vert_max) {
              /* the flip would increase the interior number of vertices over the max, cancel */
              trace("- i={}, flip cut canceled, interior nb of vertices would be too high", i);
              continue;
            }
            if (opt.val_int_max != ID_MAX && qVal[v1] >= opt.val_int_max && qVal[v2] >= opt.val_int_max) {
              /* the flip would increase the interior valence over the max, cancel */
              trace("- i={}, flip cut canceled, interior valence would be too high", i);
              continue;
            }
            if (!is_compatible(Q_edges,Q_diags,new_quad)) {
              trace("- i={}, flip cut canceled, incompatible", i);
              continue; /* compat. test */
            }
            if (i < M.loop.size()-1) {
              new_bdrLoop.erase(new_bdrLoop.begin() + i, new_bdrLoop.begin() + i + 2);
            } else {
              new_bdrLoop.erase(new_bdrLoop.begin() + i); /* last vertex of the loop */
              new_bdrLoop.erase(new_bdrLoop.begin()); /* first vertex of the loop */
            }
            n2 = M.n;
            trace("- i={}, flip cut, vPrev={}, v1={}, v2={}, vNext={}, new_quad={}", i, vPrev, v1, v2, vNext, new_quad);
          } else {
            continue;
          }
          if (new_bdrLoop.size() < 4) {
            trace("- cancel new mesh, new_bdrLoop = {}", new_bdrLoop);
            continue;
          }

          new_quad = reorder_quad(new_quad);

          QMesh M2;
          M2.n = n2;
          M2.loop = new_bdrLoop;
          M2.quads = M.quads;
          M2.quads.push_back(new_quad);

          bool add_to_meshes = true;
          bool add_to_queue = true;
          process_mesh_stats(M2, opt, add_to_meshes, add_to_queue);
          if (add_to_meshes || add_to_queue) {
            M2.loop.shrink_to_fit();
            M2.quads.shrink_to_fit();

            meshes.push_back(M2);
            meshes_add_store.push_back(add_to_meshes);
            meshes_add_queue.push_back(add_to_queue);
          }
        }
      } /* end of loop over flips */

      /* Add meshes to storage */
      trace("try to add the {} new meshes built by flip",meshes.size());
      F(i,meshes.size()) {
        if (meshes_add_store[i]) {
          store.add_quad_mesh_if_no_equivalent(meshes[i],!opt.no_nauty_check);
        }
        if (meshes_add_queue[i]) {
          bool isNew = explored.add_quad_mesh_if_no_equivalent(meshes[i],!opt.no_nauty_check);
          if (isNew) {
            queue.push(meshes[i]);
          }
        }
      }
      meshes.clear();
      meshes_add_queue.clear();
      meshes_add_store.clear();

      explored.add_quad_mesh_if_no_equivalent(M,!opt.no_nauty_check);
    } /* end of while */

    info("Found {} unique quad meshes", store.meshes.size());
    F(i,store.meshes.size()) {
      makeValenceCanonical(store.meshes[i]);
    }
    if (opt.verbose || opt.cpp) {
      F(i,store.meshes.size()) {
        const QMesh& M = store.meshes[i];
        // info("- {}/{}: {}",i+1,store.meshes.size(),string_valence_repartition(M));
        if (opt.cpp) {
          std::string str = fmt::format("{}",M.quads);
          str = ReplaceAll(str, "[", "{");
          str = ReplaceAll(str, "]", "}");
          info("  quads: {}", str);
        }
      }
    }

    if (opt.output_file.size()) {
      export_quadrangulations(store.meshes,opt.output_file,opt.cpp,opt.manual);
    }

    if (opt.visualization) {
      gmsh::initialize();
      visualize_quadrangulations(store.meshes);
      // visualize_concave_quadrangulations(store.meshes);
      gmsh::fltk::run();
      gmsh::finalize();
    }

    return true;
  }

}
