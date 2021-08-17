#ifndef BITTREE_MORTONTREE_H__
#define BITTREE_MORTONTREE_H__

#include "Bittree_prelude.h"
#include "Bittree_bitarray.h"
#include "Bittree_ref.h"

#include <iostream>

namespace BitTree {
  template<unsigned D>
  unsigned rect_coord_to_mort(const unsigned domain[D], const unsigned coord[D]);
  
  template<unsigned D>
  void rect_mort_to_coord(const unsigned domain[D], unsigned mort, unsigned coord[D]);
  
  template<unsigned D, class W>
  class MortonTree {
    unsigned levs;
    unsigned lev0_blks[D];
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
      X coord[D];
    };
  public:
    static Ref<MortonTree<D,W> > make(
      const unsigned blks[D],
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
    bool inside(unsigned lev, const X coord[D]) const;
    template<class X>
    Block<X> identify(unsigned lev, const X coord[D]) const;
    template<class X>
    Block<X> locate(X id) const;
    
    Ref<MortonTree<D,W> > refine(Ref<BitArray > delta) const;

    void bitid_list(unsigned mort_min,unsigned mort_max, int *out ) const;
    void print_if_2d(unsigned datatype=0) const;
  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;
  };

  template<class W>
  std::ostream& operator<<(std::ostream &o, const MortonTree<2,W> *x);
}
#endif
