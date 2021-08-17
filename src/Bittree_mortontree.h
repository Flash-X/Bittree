#ifndef BITTREE_MORTONTREE_H__
#define BITTREE_MORTONTREE_H__

#include "Bittree_prelude.h"
#include "Bittree_bitarray.h"
#include "Bittree_ref.h"

#include "Bittree_constants.h"

#include <iostream>

namespace BitTree {
  unsigned rect_coord_to_mort(const unsigned domain[NDIM], const unsigned coord[NDIM]);
  
  void rect_mort_to_coord(const unsigned domain[NDIM], unsigned mort, unsigned coord[NDIM]);
  
  class MortonTree {
    unsigned levs;
    unsigned lev0_blks[NDIM];
    Ref<FastBitArray > bits;
    unsigned id0; // id of first block
    struct Level {
      unsigned id1; // exclusive upper bound on block ids for this level
    } level[1];/*[levs]*/
  private:
    MortonTree() {}
  public:
    template<class X>
    struct Block {
      unsigned id;
      unsigned mort;
      unsigned level;
      bool is_parent;
      X coord[NDIM];
    };
  public:
    static Ref<MortonTree > make(
      const unsigned blks[NDIM],
      const bool includes[]=0x0
    );
    unsigned levels() const;
    unsigned blocks() const;
    unsigned leaves() const;
    unsigned top_size(unsigned dim) const;
    unsigned id_upper_bound() const;
    unsigned level_id0(unsigned lev) const;
    unsigned level_blocks(unsigned lev) const;
    void level_ids(unsigned lev, int *ids) const;
    
    bool block_is_parent(unsigned id) const;
    unsigned block_level(unsigned id) const;
    
    template<class X>
    bool inside(unsigned lev, const X coord[NDIM]) const;
    template<class X>
    Block<X> identify(unsigned lev, const X coord[NDIM]) const;
    template<class X>
    Block<X> locate(X id) const;
    
    Ref<MortonTree > refine(Ref<BitArray > delta) const;

    void bitid_list(unsigned mort_min,unsigned mort_max, int *out ) const;
    void print_if_2d(unsigned datatype=0) const;
  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;
  };

  template<class W>
  std::ostream& operator<<(std::ostream &o, const MortonTree *x);
}
#endif
