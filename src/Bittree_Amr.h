#ifndef BITTREE_AMR_H__
#define BITTREE_AMR_H__

#include "Bittree_bitarray.h"
#include "Bittree_mortontree.h"
#include "mpi.h"

namespace bittree {
  
class BittreeAmr  {
  public:
    BittreeAmr(const int top[], const int includes[]);

    std::shared_ptr<MortonTree> getTree(bool updated=false);

    // Get refinement info
    unsigned delta_count() const;
    bool check_refine_bit(unsigned bitid) const;

    // Refinement functions
    void refine_init();
    void refine_mark(unsigned bitid, bool value);
    void refine_reduce(MPI_Comm comm);
    void refine_reduce_and(MPI_Comm comm);
    void refine_update();
    void refine_apply();

    // Other functions
    std::string slice_to_string(unsigned datatype, unsigned slice=0) const;

  private:
    std::shared_ptr<MortonTree> tree_;            //!<Actual Bittree
    std::shared_ptr<MortonTree> tree_updated_;    //!<Updated Bittree, before refinement is applied
    std::shared_ptr<BitArray> refine_delta_;      //!<(De)refinement flags for blocks
    bool is_reduced_;  //!<Flag to track whether refine_delta is up to date across processors
    bool is_updated_;  //!<Flag to track whether tree_updated matches latest refine_delta
    bool in_refine_;   //!<If in_refine=false, tree_updated and refine_delta should not exist
  };

}

#endif
