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
  
  std::shared_ptr<MortonTree > MortonTree::make(
      const unsigned size[NDIM],
      const bool includes[]
    ) {
    std::shared_ptr<MortonTree> it = std::make_shared<MortonTree>();
    
    unsigned blkpop = 1;
    for(unsigned d=0; d < NDIM; d++) {
      it->lev0_blks[d] = size[d];
      blkpop *= size[d];
    }
    it->levs = 1;
    // the first blkpop bits of our bitarray are 'inclusion' bits
    // after that come the actual block bits for all but the last level, which has no bits
    it->id0 = blkpop;
    unsigned lev0_id1 = it->id0;
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
    it->bits = bldr.finish();
    it->level.push_back(LevelStruct{.id1 = lev0_id1});
    return it;
  }

  unsigned MortonTree::levels() const {
    return this->levs;
  }

  unsigned MortonTree::blocks() const {
    return this->level[this->levs-1].id1 - this->id0;
  }

  unsigned MortonTree::leaves() const {
    unsigned pars = bits->count( this->id0, this->level[this->levs-1].id1) ;
    return this->blocks() - pars ;
  }

  unsigned MortonTree::top_size(unsigned dim) const {
    return this->lev0_blks[dim];
  }
  
  unsigned MortonTree::id_upper_bound() const {
    return this->level[this->levs-1].id1;
  }

  unsigned MortonTree::level_id0(unsigned lev) const {
    return lev == 0 ? this->id0 : this->level[lev-1].id1;
  }

  unsigned MortonTree::level_blocks(unsigned lev) const {
    return this->level[lev].id1 - (lev == 0 ? this->id0 : this->level[lev-1].id1);
  }

  void MortonTree::level_ids(unsigned lev, int* ids) const {
    ids[1] = static_cast<int>( this->level[lev].id1 );
    ids[0] = static_cast<int>( lev == 0 ? this->id0 : this->level[lev-1].id1 ); 
    return;
  }

  bool MortonTree::block_is_parent(unsigned id) const {
    if(levs>1) return id < level[levs-2].id1 && bits->get(id);
    else return false;
  }

  unsigned MortonTree::block_level(unsigned id) const {
    unsigned lev = 0;
    while(level[lev].id1 <= id)
      lev += 1;
    return lev;
  }

  bool MortonTree::inside(unsigned lev, const unsigned x[NDIM]) const {
    unsigned x0[NDIM];
    for(unsigned d=0; d < NDIM; d++) {
      x0[d] = x[d] >> lev;
      if(x0[d] >= lev0_blks[d])
        return false;
    }
    return bits->get(rect_coord_to_mort(lev0_blks, x0));
  }
  
  /** Identifies morton number of a block corresponding to given coords
   *  on the current tree.  */
  MortonTree::Block
  MortonTree::identify(unsigned lev, const unsigned x[NDIM]) const {
    const unsigned levs = this->levs;
    const unsigned id0 = this->id0;
    const std::shared_ptr<BitArray> bits_a = bits->bit_array(); // use this for bit access
    Block ans;
    unsigned ix; // index of current block in current level
    { // top level=0
      unsigned x0[NDIM];
      for(unsigned d=0; d < NDIM; d++) {
        x0[d] = unsigned(x[d] >> lev); // coarsen x to top level
        DBG_ASSERT(x0[d] < lev0_blks[d]);
        ans.coord[d] = unsigned(x0[d]);
      }
      ix = rect_coord_to_mort(lev0_blks, x0);
      DBG_ASSERT(bits_a->get(ix));
      ix = bits->count(0, ix); // discount excluded blocks
    }
    ans.mort = 0;
    // bisection iteration
    for(unsigned a_lev=0; a_lev < levs; a_lev++) {
      unsigned a_id0 = a_lev==0 ? id0 : level[a_lev-1].id1;
      unsigned a_id = a_id0 + ix;
      ans.mort += ix;
      unsigned inside = 0u;
      bool is_par = a_lev+1u < levs && bits_a->get(a_id);
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
      unsigned parbef = this->parents_before(a_lev, ix);
      ix = (parbef<<NDIM) + inside;
    }
    return ans;
  }

  MortonTree::Block
  MortonTree::locate(unsigned id) const {
    const unsigned id0 = this->id0;
    Block ans;
    ans.id = id;
    ans.mort = 0;
    ans.is_parent = this->block_is_parent(id);
    unsigned lev = 0;
    while(level[lev].id1 <= id)
      lev += 1;
    ans.level = lev;
    for(unsigned d=0; d < NDIM; d++)
      ans.coord[d] = unsigned(0u);
    // index on this level
    unsigned ix = id - (lev == 0 ? id0 : level[lev-1].id1);
    { // count children of all preceeding parents in morton index
      unsigned down = ix;
      for(unsigned lev1=lev; lev1 < levs; lev1++) {
        down = this->parents_before(lev1, down) << NDIM;
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
      ix = this->parent_find(lev-1, ix>>NDIM) - (lev-1==0 ? id0 : level[lev-2].id1);
      lev -= 1;
    }
    ans.mort += ix;
    { // top level=0
      unsigned x0[NDIM];
      ix = bits->find(0, ix); // account for excluded blocks
      rect_mort_to_coord(lev0_blks, ix, x0);
      for(unsigned d=0; d < NDIM; d++)
        ans.coord[d] += unsigned(x0[d]) << ans.level;
    }
    return ans;
  }

  std::shared_ptr<MortonTree > MortonTree::refine(
      std::shared_ptr<BitArray > delta_
    ) const {
    
    const unsigned a_levs = this->levs;
    const std::shared_ptr<BitArray> a_bits = this->bits->bit_array();
    const unsigned id0 = this->id0;
    
    // count the new number of levels, blocks, and bits
    unsigned b_id1 = this->level[0].id1;
    unsigned b_bitlen = id0;
    unsigned b_levs = 1;
    for(unsigned lev=0; lev < a_levs; lev++) {
      unsigned lev_id0 = lev == 0 ? id0 : this->level[lev-1].id1;
      unsigned lev_id1 = this->level[lev].id1;
      unsigned b_pars = BitArray::count_xor(*a_bits, *delta_, lev_id0, lev_id1);
      if(b_pars != 0) b_bitlen = b_id1;
      b_id1 += b_pars << NDIM;
      if(b_pars == 0) break;
      b_levs += 1;
    }
    
    // new bit tree
    std::shared_ptr<MortonTree> b_tree = std::make_shared<MortonTree>();
    {
      b_tree->levs = b_levs;
      b_tree->id0 = id0;
      b_tree->level.resize(b_levs);
      for(unsigned d=0; d < NDIM; d++)
        b_tree->lev0_blks[d] = this->lev0_blks[d];
      // still must initialize b_tree->bits
    }
    
    // apply delta and insert/remove blocks
    typename BitArray::Reader a_r(a_bits), del_r(delta_, id0);
    typename FastBitArray::Builder b_w(b_bitlen);
    
    // copy inclusion bits
    while(a_r.index() < id0)
      b_w.write<1>(a_r.read<1>());
    
    // do level 0...
    b_tree->level[0].id1 = this->level[0].id1;
    while(b_w.index() < b_bitlen && a_r.index() < this->level[0].id1) {
      // apply delta
      b_w.write<1>(a_r.read<1>() ^ del_r.read<1>());
    }
    
    // readers of previous level
    typename BitArray::Reader a_rp(a_bits, id0), del_rp(delta_, id0);

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
      if(a_rp.index() == this->level[lev-1].id1 || b_w.index() == b_bitlen) {
        b_tree->level[lev].id1 = b_w.index();
        lev += 1;
      }
    }
    
    b_tree->level[b_levs-1].id1 = b_id1;
    b_tree->bits = b_w.finish();

    return b_tree;
  }

  unsigned MortonTree::parents_before(unsigned lev, unsigned ix) const {
    if(lev >= levs-1) return 0;
    unsigned id0 = lev == 0 ? this->id0 : level[lev-1].id1;
    return bits->count(id0, id0 + ix);
  }

  unsigned MortonTree::parent_find(unsigned lev, unsigned par_ix) const {
    unsigned id0 = lev == 0 ? this->id0 : level[lev-1].id1;
    return bits->find(id0, par_ix);
  }

  void MortonTree::bitid_list(unsigned mort_min, unsigned mort_max, int *out ) const {
    bool is_par; 
    unsigned ix = this->id0;           //current scan index
    unsigned lev = 0;          //current scanning level
    bool childrenDone[this->levs];
    unsigned pos[this->levs]; //location on each level (increases monotonically)
    unsigned mort = 0;
#ifdef FLASH_DEBUG_BITTREE
    DBG_ASSERT(mort_max <= this->blocks());
    DBG_ASSERT(mort_min <= mort_max);
#endif

    for (unsigned i=0; i<levs; i++){
      childrenDone[i] = false;
      pos[i] = 0;
    }

    //each iteration either advances ix by one, or goes up/down one level
    while(mort < mort_max) {
      is_par = this->block_is_parent(ix);

      //if scanning a parent and children have not been scanned, move down a level
      if(is_par && !childrenDone[lev]) {
        ix = this->level[lev].id1 + ((1u<<NDIM) * this->parents_before(lev,pos[lev]));
        childrenDone[lev+1]=false;
       
#ifndef ALT_MORTON_ORDER
        if(mort<mort_max && mort>=mort_min) out[mort-mort_min] = int(pos[lev] + this->level_id0(lev)) ;
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
          if(mort<mort_max && mort>=mort_min) out[mort-mort_min] = int(pos[lev-1] + this->level_id0(lev-1)) ;
          mort++;
        }
#endif

        //if last child
        if (lev>0 && (((pos[lev]+1) % (1u<<NDIM)) == 0) ) {
          pos[lev]++;
          childrenDone[lev-1] = true;
          ix = pos[lev-1] + this->level_id0(lev-1);
          lev--;
        }
        //if last block on top level
        else if (lev==0 && (pos[lev]+1)==this->level[0].id1 ){
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
    
    unsigned levs = this->levels();
    for(unsigned lev=0; lev < levs; lev++) {
      cout << "lev=" << lev << '\n';
      
      unsigned ij[2];
      for(ij[1]=0; ij[1] < this->top_size(1)<<lev; ij[1]++) {
        for(ij[0]=0; ij[0] < this->top_size(0)<<lev; ij[0]++) {
          if(this->inside(lev, ij)) {
            MortonTree::Block b0 = this->identify(lev, ij);
            MortonTree::Block b1 = this->locate(b0.id);
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
                cout << std::setw(4) << this->block_is_parent(b0.id); //print block parentage
                break;
            }
          }
          else
            cout << std::setw(4) << ' ';
        }
        cout << '\n';
      }
    }
  }

  std::ostream& operator<<(std::ostream &o, const MortonTree *x) {
    using namespace std;
    
    unsigned levs = x->levels();
    for(unsigned lev=0; lev < levs; lev++) {
      cout << "lev=" << lev << '\n';
      unsigned ij[2];
      for(ij[1]=0; ij[1] < x->top_size(1)<<lev; ij[1]++) {
        for(ij[0]=0; ij[0] < x->top_size(0)<<lev; ij[0]++) {
          if(x->inside(lev, ij)) {
            MortonTree::Block b0 = x->identify(lev, ij);
            MortonTree::Block b1 = x->locate(b0.id);
            DBG_ASSERT(b0.id == b1.id);
            DBG_ASSERT(b0.level == b1.level);
            DBG_ASSERT(b0.mort == b1.mort);
            for(unsigned i=0; i<NDIM; ++i) {
              DBG_ASSERT(b0.coord[i] == b1.coord[i]);
            }
            o << std::setw(4) << b0.mort;
          }
          else
            o << std::setw(4) << ' ';
        }
        o << '\n';
      }
    }
    return o;
  }
}
