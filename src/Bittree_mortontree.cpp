#include "Bittree_mortontree.h"

#include "Bittree_bitarray.h"
#include "Bittree_bits.h"

#include <algorithm>
#include <limits>
#include <iomanip>
#include <iostream>

namespace bittree {
  unsigned rect_coord_to_mort(const unsigned domain[NDIM], const unsigned coord[NDIM]) {
    unsigned x[NDIM], box[NDIM];
    for(unsigned d=0; d < NDIM; d++) {
      x[d] = coord[d];
      box[d] = domain[d];
    }
    unsigned mort = 0u;
    // bisect box until it's just one element
    while(true) {
      // find dim that can fit biggest pow2 strictly inside box
      unsigned max_pow2 = 0u;
      unsigned max_d;
      for(unsigned d=0; d < NDIM; d++) {
        unsigned p2 = glb_pow2(box[d]-1u);
        if(p2 >= max_pow2) {
          max_pow2 = p2;
          max_d = d;
        }
      }
      if(max_pow2 == 0u)
        return mort; // the box is just one, we're done
      if(x[max_d] < max_pow2)
        box[max_d] = max_pow2;
      else {
        unsigned pop = 1u;
        for(unsigned d=0; d < NDIM; d++)
          pop *= d == max_d ? max_pow2 : box[d];
        mort += pop;
        x[max_d] -= max_pow2;
        box[max_d] -= max_pow2;
      }
    }
  }
  
  void rect_mort_to_coord(const unsigned domain[NDIM], unsigned mort, unsigned coord[NDIM]) {
    unsigned box[NDIM];
    for(unsigned d=0; d < NDIM; d++) {
      coord[d] = 0u;
      box[d] = domain[d];
    }
    // bisect box until it's just one element
    while(true) {
      // find dim that can fit biggest pow2 strictly inside box
      unsigned max_pow2 = 0u;
      unsigned max_d;
      for(unsigned d=0; d < NDIM; d++) {
        unsigned p2 = glb_pow2(box[d]-1u);
        if(p2 >= max_pow2) {
          max_pow2 = p2;
          max_d = d;
        }
      }
      if(max_pow2 == 0u)
        return; // the box is just one, we're done
      unsigned pop = 1u; // left sub-box population
      for(unsigned d=0; d < NDIM; d++)
        pop *= d == max_d ? max_pow2 : box[d];
      if(mort < pop)
        box[max_d] = max_pow2;
      else {
        mort -= pop;
        coord[max_d] += max_pow2;
        box[max_d] -= max_pow2;
      }
    }
  }
  
  MortonTree::MortonTree(const unsigned size[NDIM], const bool includes[]) {
    
    unsigned blkpop = 1;
    for(unsigned d=0; d < NDIM; d++) {
      lev0_blks_[d] = size[d];
      blkpop *= size[d];
    }
    levs_ = 1;
    // the first blkpop bits of our bitarray are 'inclusion' bits
    // after that come the actual block bits for all but the last level, which has no bits
    id0_ = blkpop;
    unsigned lev0_id1 = id0_;
    typename FastBitArray::Builder bldr(blkpop);
    // generate inclusion bits
    for(unsigned mort=0; mort < blkpop; mort++) {
      unsigned x[NDIM];
      rect_mort_to_coord(size, mort, x);
      unsigned ix = 0;
      for(unsigned d=NDIM; d--;)
        ix = size[d]*ix + x[d];
      unsigned include = includes[ix] ? 1 : 0;
      lev0_id1 += include;
      bldr.write<1>(include);
    }
    // ok bitarray done.  since there's only one level, we dont store any block bits
    bits_ = bldr.finish();
    level_.push_back(LevelStruct{.id1 = lev0_id1});
  }

  unsigned MortonTree::levels() const {
    return levs_;
  }

  unsigned MortonTree::blocks() const {
    return level_[levs_-1].id1 - id0_;
  }

  unsigned MortonTree::leaves() const {
    unsigned pars = bits_->count( id0_, level_[levs_-1].id1) ;
    return blocks() - pars ;
  }

  unsigned MortonTree::top_size(unsigned dim) const {
    return lev0_blks_[dim];
  }
  
  unsigned MortonTree::id_upper_bound() const {
    return level_[levs_-1].id1;
  }

  /**
    * \todo Inline this function?
    */
  unsigned MortonTree::level_id0(unsigned lev) const {
    return lev == 0 ? id0_ : level_[lev-1].id1;
  }

  unsigned MortonTree::level_blocks(unsigned lev) const {
    return level_[lev].id1 - (lev == 0 ? id0_ : level_[lev-1].id1);
  }

  /**
    * \todo Adjust this interface - this is more of a Fortran layer thing
    */
  void MortonTree::level_ids(unsigned lev, int* ids) const {
    ids[0] = static_cast<int>( level_id0(lev) );
    ids[1] = static_cast<int>( level_[lev].id1 );
    return;
  }

  bool MortonTree::block_is_parent(unsigned id) const {
    if(levs_>1) return id < level_[levs_-2].id1 && bits_->get(id);
    else return false;
  }

  unsigned MortonTree::block_level(unsigned id) const {
    unsigned lev = 0;
    while(level_[lev].id1 <= id)
      lev += 1;
    return lev;
  }

  bool MortonTree::inside(unsigned lev, const unsigned x[NDIM]) const {
    unsigned x0[NDIM];
    for(unsigned d=0; d < NDIM; d++) {
      x0[d] = x[d] >> lev;
      if(x0[d] >= lev0_blks_[d])
        return false;
    }
    return bits_->get(rect_coord_to_mort(lev0_blks_, x0));
  }
  
  /** Identifies morton number of a block corresponding to given coords
   *  on the current tree.  */
  MortonTree::Block
  MortonTree::identify(unsigned lev, const unsigned x[NDIM]) const {
    const std::shared_ptr<BitArray> bits_a = bits_->bit_array(); // use this for bit access
    Block ans;
    unsigned ix; // index of current block in current level
    { // top level=0
      unsigned x0[NDIM];
      for(unsigned d=0; d < NDIM; d++) {
        x0[d] = unsigned(x[d] >> lev); // coarsen x to top level
        DBG_ASSERT(x0[d] < lev0_blks_[d]);
        ans.coord[d] = unsigned(x0[d]);
      }
      ix = rect_coord_to_mort(lev0_blks_, x0);
      DBG_ASSERT(bits_a->get(ix));
      ix = bits_->count(0, ix); // discount excluded blocks
    }
    ans.mort = 0;
    // bisection iteration
    for(unsigned a_lev=0; a_lev < levs_; a_lev++) {
      unsigned a_id0 = a_lev==0 ? id0_ : level_[a_lev-1].id1;
      unsigned a_id = a_id0 + ix;
      ans.mort += ix;
      unsigned inside = 0u;
      bool is_par = a_lev+1u < levs_ && bits_a->get(a_id);
      if(is_par && a_lev < lev) { // if we're bisecting further
#ifndef ALT_MORTON_ORDER
        ans.mort += 1;
#endif
        for(unsigned d=0; d < NDIM; d++) {
          unsigned xd = x[d] >> (lev-a_lev-1u);
          ans.coord[d] <<= 1;
          if(xd >= ans.coord[d]+1u) {
            ans.coord[d] += 1u;
            inside += 1u << d;
#ifdef ALT_MORTON_ORDER
            ans.mort += d == NDIM-1 ? 1 : 0;
#endif
          }
        }
      }
      else if(a_lev <= lev) { // we have the result block
        lev = a_lev; // stop this from running again
        ans.id = a_id;
        ans.level = a_lev;
        ans.is_parent = is_par;
#ifdef ALT_MORTON_ORDER
        if(is_par) inside = 1u<<(NDIM-1); //include first half of children
#endif
      }
      unsigned parbef = parents_before(a_lev, ix);
      ix = (parbef<<NDIM) + inside;
    }
    return ans;
  }

  MortonTree::Block
  MortonTree::locate(unsigned id) const {
    Block ans;
    ans.id = id;
    ans.mort = 0;
    ans.is_parent = block_is_parent(id);
    unsigned lev = 0;
    while(level_[lev].id1 <= id)
      lev += 1;
    ans.level = lev;
    for(unsigned d=0; d < NDIM; d++)
      ans.coord[d] = unsigned(0u);
    // index on this level
    unsigned ix = id - (lev == 0 ? id0_ : level_[lev-1].id1);
    { // count children of all preceeding parents in morton index
      unsigned down = ix;
      for(unsigned lev1=lev; lev1 < levs_; lev1++) {
        down = parents_before(lev1, down) << NDIM;
#ifdef ALT_MORTON_ORDER
        if(lev1 == lev && ans.is_parent)
          down += 1u<<(NDIM-1);
#endif
        ans.mort += down;
      }
    }
    // walk up the levels
    while(0 < lev) {
      for(unsigned d=0; d < NDIM; d++)
        ans.coord[d] += unsigned(ix>>d & 1u) << (ans.level-lev);
#ifdef ALT_MORTON_ORDER
      ans.mort += ix + (ix>>(NDIM-1) & 1u);
#else
      ans.mort += ix + 1;
#endif
      ix = parent_find(lev-1, ix>>NDIM) - (lev-1==0 ? id0_ : level_[lev-2].id1);
      lev -= 1;
    }
    ans.mort += ix;
    { // top level=0
      unsigned x0[NDIM];
      ix = bits_->find(0, ix); // account for excluded blocks
      rect_mort_to_coord(lev0_blks_, ix, x0);
      for(unsigned d=0; d < NDIM; d++)
        ans.coord[d] += unsigned(x0[d]) << ans.level;
    }
    return ans;
  }

  std::shared_ptr<MortonTree > MortonTree::refine(std::shared_ptr<const BitArray> delta) const {
    
    const unsigned a_levs = levs_;
    const std::shared_ptr<BitArray> a_bits = bits_->bit_array();
    
    // count the new number of levels, blocks, and bits
    unsigned b_id1 = level_[0].id1;
    unsigned b_bitlen = id0_;
    unsigned b_levs = 1;
    for(unsigned lev=0; lev < a_levs; lev++) {
      unsigned lev_id0 = lev == 0 ? id0_ : level_[lev-1].id1;
      unsigned lev_id1 = level_[lev].id1;
      unsigned b_pars = BitArray::count_xor(*a_bits, *delta, lev_id0, lev_id1);
      if(b_pars != 0) b_bitlen = b_id1;
      b_id1 += b_pars << NDIM;
      if(b_pars == 0) break;
      b_levs += 1;
    }
    
    // new bit tree
    std::shared_ptr<MortonTree> b_tree = std::make_shared<MortonTree>();
    {
      b_tree->levs_ = b_levs;
      b_tree->id0_ = id0_;
      b_tree->level_.resize(b_levs);
      for(unsigned d=0; d < NDIM; d++)
        b_tree->lev0_blks_[d] = lev0_blks_[d];
      // still must initialize b_tree->bits
    }
    
    // apply delta and insert/remove blocks
    typename BitArray::Reader a_r(a_bits), del_r(delta, id0_);
    typename FastBitArray::Builder b_w(b_bitlen);
    
    // copy inclusion bits
    while(a_r.index() < id0_)
      b_w.write<1>(a_r.read<1>());
    
    // do level 0...
    b_tree->level_[0].id1 = level_[0].id1;
    while(b_w.index() < b_bitlen && a_r.index() < level_[0].id1) {
      // apply delta
      b_w.write<1>(a_r.read<1>() ^ del_r.read<1>());
    }
    
    // readers of previous level
    typename BitArray::Reader a_rp(a_bits, id0_), del_rp(delta, id0_);

    // do remaining levels
    unsigned lev = 1;
    while(b_w.index() < b_bitlen) {
      if(a_rp.read<1>()) { // it was a parent
        bool still_a_parent = !del_rp.read<1>();
        // read kids, apply delta
        BitArray::WType b_kids = a_r.read<(1<<NDIM)>() ^ del_r.read<(1<<NDIM)>();
        // if it became a leaf then we just dont write out the kids
        if(still_a_parent)
          b_w.write<(1<<NDIM)>(b_kids);
      }
      else { // it was a leaf
        if(del_rp.read<1>()) // and it became a parent!
          b_w.write<(1<<NDIM)>(0);
      }
      if(a_rp.index() == level_[lev-1].id1 || b_w.index() == b_bitlen) {
        b_tree->level_[lev].id1 = b_w.index();
        lev += 1;
      }
    }
    
    b_tree->level_[b_levs-1].id1 = b_id1;
    b_tree->bits_ = b_w.finish();

    return b_tree;
  }

  unsigned MortonTree::parents_before(unsigned lev, unsigned ix) const {
    if(lev >= levs_-1) return 0;
    return bits_->count(level_id0(lev), level_id0(lev) + ix);
  }

  unsigned MortonTree::parent_find(unsigned lev, unsigned par_ix) const {
    return bits_->find(level_id0(lev), par_ix);
  }

  void MortonTree::bitid_list(unsigned mort_min, unsigned mort_max, int *out ) const {
    bool is_par; 
    unsigned ix = id0_;           //current scan index
    unsigned lev = 0;          //current scanning level
    bool childrenDone[levs_];
    unsigned pos[levs_]; //location on each level (increases monotonically)
    unsigned mort = 0;
#ifdef FLASH_DEBUG_BITTREE
    DBG_ASSERT(mort_max <= blocks());
    DBG_ASSERT(mort_min <= mort_max);
#endif

    for (unsigned i=0; i<levs_; i++){
      childrenDone[i] = false;
      pos[i] = 0;
    }

    //each iteration either advances ix by one, or goes up/down one level
    while(mort < mort_max) {
      is_par = block_is_parent(ix);

      //if scanning a parent and children have not been scanned, move down a level
      if(is_par && !childrenDone[lev]) {
        ix = level_[lev].id1 + ((1u<<NDIM) * parents_before(lev,pos[lev]));
        childrenDone[lev+1]=false;
       
#ifndef ALT_MORTON_ORDER
        if(mort<mort_max && mort>=mort_min) out[mort-mort_min] = int(pos[lev] + level_id0(lev)) ;
        mort++;
#endif

        lev++;
      }
      else {
        //if leaf, store its bitid
        if (!is_par) {
          if(mort<mort_max && mort>=mort_min) out[mort-mort_min] = int(ix);
          mort++;
        }

#ifdef ALT_MORTON_ORDER
        //if middle child, store parent's bitid
        if (lev>0 && (((pos[lev]+1) % (1u<<NDIM)) == (1u<<(NDIM-1))) ){
          if(mort<mort_max && mort>=mort_min) out[mort-mort_min] = int(pos[lev-1] + level_id0(lev-1)) ;
          mort++;
        }
#endif

        //if last child
        if (lev>0 && (((pos[lev]+1) % (1u<<NDIM)) == 0) ) {
          pos[lev]++;
          childrenDone[lev-1] = true;
          ix = pos[lev-1] + level_id0(lev-1);
          lev--;
        }
        //if last block on top level
        else if (lev==0 && (pos[lev]+1)==level_[0].id1 ){
          return;
        }
        //else
        else {
          pos[lev]++;
          ix++;
          childrenDone[lev] = false;
        }
      }
    }
  }

  void MortonTree::print_if_2d(unsigned datatype) const {
    DBG_ASSERT(NDIM==2);
    using namespace std;
    
    unsigned levs = levels();
    for(unsigned lev=0; lev < levs; lev++) {
      cout << "lev=" << lev << '\n';
      
      unsigned ij[2];
      for(ij[1]=0; ij[1] < top_size(1)<<lev; ij[1]++) {
        for(ij[0]=0; ij[0] < top_size(0)<<lev; ij[0]++) {
          if(inside(lev, ij)) {
            MortonTree::Block b0 = identify(lev, ij);
            MortonTree::Block b1 = locate(b0.id);
            DBG_ASSERT(b0.id == b1.id);
            DBG_ASSERT(b0.level == b1.level);
            DBG_ASSERT(b0.mort == b1.mort);
            DBG_ASSERT(b0.coord[0] == b1.coord[0] && b0.coord[1] == b1.coord[1]);
            switch(datatype) {
              case 0:
                cout << std::setw(4) << b0.id;         //print bittree id number
                break;
              case 1: 
                cout << std::setw(4) << (b0.mort+1);   //print 1-based morton number
                break;
              case 2: 
                cout << std::setw(4) << block_is_parent(b0.id); //print block parentage
                break;
            }
          }
          else
            cout << std::setw(4) << ' ';
        }
        cout << '\n';
      }
    }
    cout.flush();
  }

}
