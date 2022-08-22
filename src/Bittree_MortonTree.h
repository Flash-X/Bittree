/*
   Copyright 2022 UChicago Argonne, LLC and contributors

   Licensed under the Apache License, Version 2.0 (the "License"); 
   you may not use this file except in compliance with the License. 
    
 
   Unless required by applicable law or agreed to in writing, software 
   distributed under the License is distributed on an "AS IS" BASIS, 
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
   See the License for the specific language governing permissions and 
   limitations under the License.
*/
#ifndef BITTREE_MORTONTREE_H__
#define BITTREE_MORTONTREE_H__

#include "Bittree_BitArray.h"
#include "Bittree_constants.h"

namespace bittree {
  unsigned rect_coord_to_mort(const unsigned domain[BTDIM], const unsigned coord[BTDIM]);
  void rect_mort_to_coord(const unsigned domain[BTDIM], unsigned mort, unsigned coord[BTDIM]);

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
      unsigned coord[BTDIM];
    };

    struct LevelStruct {
      unsigned id1; // exclusive upper bound on block ids for this level
    };


  public:
    MortonTree() {}
    MortonTree(const int blks[BTDIM], const int includes[]);
    ~MortonTree() = default;

    // Getters
    unsigned levels() const;
    unsigned blocks() const;
    unsigned leaves() const;
    unsigned top_size(unsigned dim) const;
    unsigned id_upper_bound() const;
    unsigned level_id0(unsigned lev) const;
    unsigned level_id1(unsigned lev) const;
    unsigned level_blocks(unsigned lev) const;

    // Other member functions
    unsigned getParentId(unsigned id) const;
    bool block_is_parent(unsigned id) const;
    unsigned block_level(unsigned id) const;
    Block locate(unsigned id) const;
    bool inside(unsigned lev, const unsigned coord[BTDIM]) const;
    Block identify(unsigned lev, const unsigned coord[BTDIM]) const;

    std::shared_ptr<MortonTree> refine(std::shared_ptr<const BitArray> delta) const;
    void bitid_list(unsigned mort_min,unsigned mort_max, int *out ) const;

    std::string print_slice(unsigned datatype, unsigned slice=0) const;

  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;

  public:
    std::shared_ptr<FastBitArray> bits_;   //!< Data
  private:
    // Member variables
    unsigned levs_;                        //!< Current number of levels
    unsigned lev0_blks_[BTDIM];             //!< Number of top level blocks
    unsigned id0_;                         //!< id of first block
    std::vector<LevelStruct> level_;       //!< Upper bound on ids for each level
  };

}
#endif
