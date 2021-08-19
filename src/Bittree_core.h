#ifndef BITTREE_CORE_H__
#define BITTREE_CORE_H__

#include "Bittree_bitarray.h"
#include "Bittree_mortontree.h"
#include "mpi.h"
#include <iostream>

using namespace bittree;
  
  class TheTree  {
  public:
    TheTree(const unsigned top[], const bool includes[]);

    unsigned level_count(bool updated);
    unsigned block_count(bool updated);
    unsigned leaf_count(bool updated);
    unsigned delta_count();
    bool check_refine_bit(unsigned bitid);
    bool is_parent(bool updated, unsigned bitid);
    void identify(bool updated, int *lev, int *ijk, int *mort, int *bitid);
    void locate(bool updated, unsigned bitid, int *lev, int *ijk, int *mort);
    void get_id0(bool updated, int *out);
    void get_level_id_limits(bool updated, unsigned lev, int *ids);
    void get_bitid_list(bool updated, unsigned mort_min, unsigned mort_max, int *out);
    void refine_init();
    void refine_mark(unsigned bitid, bool value);
    void refine_reduce(MPI_Comm comm);
    void refine_reduce_and(MPI_Comm comm);
    void refine_update();
    void refine_apply();
    void print_2d(unsigned datatype=0);

  private:
    std::shared_ptr<MortonTree> tree_;            //!<Actual Bittree
    std::shared_ptr<MortonTree> tree_updated_;    //!<Updated Bittree, before refinement is applied
    std::shared_ptr<BitArray> refine_delta_;      //!<(De)refinement flags for blocks
    bool is_reduced_;  //!<Flag to track whether refine_delta is up to date across processors
    bool is_updated_;  //!<Flag to track whether tree_updated matches latest refine_delta
    bool in_refine_;   //!<If in_refine=false, tree_updated and refine_delta should not exist
  };
  

#endif
