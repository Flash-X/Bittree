#ifndef BITTREE_MORTONTREE_H__
#define BITTREE_MORTONTREE_H__

#include "Bittree_prelude.h"
#include "Bittree_bitarray.h"

#include "Bittree_constants.h"

#include <iostream>

namespace bittree {
  unsigned rect_coord_to_mort(const unsigned domain[NDIM], const unsigned coord[NDIM]);
  void rect_mort_to_coord(const unsigned domain[NDIM], unsigned mort, unsigned coord[NDIM]);

  /** MortonTree is a bitmap with metadata to make it a properly defined mesh.
   *
   */
  class MortonTree {

  public:
    struct Block {
      unsigned id;
      unsigned mort;
      unsigned level;
      bool is_parent;
      unsigned coord[NDIM];
    };

    struct LevelStruct {
      unsigned id1; // exclusive upper bound on block ids for this level
    };


  public:
    MortonTree() {}
    MortonTree(const unsigned blks[NDIM], const bool includes[]=0x0);
    ~MortonTree() = default;

    // Getters
    unsigned levels() const;
    unsigned blocks() const;
    unsigned leaves() const;
    unsigned top_size(unsigned dim) const;
    unsigned id_upper_bound() const;
    unsigned level_id0(unsigned lev) const;
    unsigned level_blocks(unsigned lev) const;
    void level_ids(unsigned lev, int *ids) const;

    // Other member functions
    bool block_is_parent(unsigned id) const;
    unsigned block_level(unsigned id) const;
    Block locate(unsigned id) const;
    bool inside(unsigned lev, const unsigned coord[NDIM]) const;
    Block identify(unsigned lev, const unsigned coord[NDIM]) const;

    std::shared_ptr<MortonTree> refine(std::shared_ptr<const BitArray> delta) const;
    void bitid_list(unsigned mort_min,unsigned mort_max, int *out ) const;

    void print_if_2d(unsigned datatype=0) const;

  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;

  private:
    // Member variables
    unsigned levs_;                        //!< Current number of levels
    unsigned lev0_blks_[NDIM];             //!< Number of top level blocks
    std::shared_ptr<FastBitArray> bits_;   //!< Data
    unsigned id0_;                         //!< id of first block
    std::vector<LevelStruct> level_;       //!< Upper bound on ids for each level
  };

}
#endif
