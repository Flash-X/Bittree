#ifndef BITTREE_MORTONTREE_H__
#define BITTREE_MORTONTREE_H__

#include "Bittree_prelude.h"
#include "Bittree_bitarray.h"

#include "Bittree_constants.h"

#include <iostream>

namespace bittree {
  struct LevelStruct {
    unsigned id1; // exclusive upper bound on block ids for this level
  };

  unsigned rect_coord_to_mort(const unsigned domain[NDIM], const unsigned coord[NDIM]);
  
  void rect_mort_to_coord(const unsigned domain[NDIM], unsigned mort, unsigned coord[NDIM]);
  
  class MortonTree {
    unsigned levs;
    unsigned lev0_blks[NDIM];
    std::shared_ptr<FastBitArray > bits;
    unsigned id0; // id of first block
    std::vector<LevelStruct> level;
  public:
    MortonTree() {}
    ~MortonTree() {std::cout << "Calling MT destructor" << std::endl;}
  public:
    struct Block {
      unsigned id;
      unsigned mort;
      unsigned level;
      bool is_parent;
      unsigned coord[NDIM];
    };
  public:
    static std::shared_ptr<MortonTree > make(
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
    
    bool inside(unsigned lev, const unsigned coord[NDIM]) const;
    Block identify(unsigned lev, const unsigned coord[NDIM]) const;
    Block locate(unsigned id) const;
    
    std::shared_ptr<MortonTree > refine(std::shared_ptr<BitArray > delta) const;

    void bitid_list(unsigned mort_min,unsigned mort_max, int *out ) const;
    void print_if_2d(unsigned datatype=0) const;
  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;
  };

  std::ostream& operator<<(std::ostream &o, const MortonTree *x);
}
#endif
