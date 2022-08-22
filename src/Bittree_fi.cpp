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
#include "Bittree_fi.h"

#include <iostream>

/** Checks if the_tree has been created */
extern "C" bool bittree_initialized() {
  return !!the_tree;
}

/** Essentially a wrapper for TheTree's constructor */
extern "C" void bittree_init(
    int topsize[],  // in
    int includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
  ) {
 
  the_tree = std::make_shared<BittreeAmr>(topsize,includes);
}

/** Wrapper function for block_count */
extern "C" void bittree_level_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    *count = static_cast<int>(tree->levels());
  }
}


/** Wrapper function for block_count */
extern "C" void bittree_block_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    *count = static_cast<int>(tree->blocks());
  }
}


/** Wrapper function for leaf_count */
extern "C" void bittree_leaf_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    *count = static_cast<int>(tree->leaves());
  }
}

/** Wrapper function for delta_count */
extern "C" void bittree_delta_count(
    int *count         //out
  ) {
  if(!!the_tree)
    *count = static_cast<int>(the_tree->delta_count());
}


/** Wrapper function for check_refine_bit */
extern "C" void bittree_check_refine_bit(
    const int *bitid,   //in
    bool *bit_check     //out
  ) {
  unsigned bitid_u  = static_cast<unsigned>(*bitid);
  if(!!the_tree)
    *bit_check = the_tree->check_refine_bit(bitid_u);
}

/** Wrapper function for is_parent */
extern "C" void bittree_is_parent(
    bool *updated,      //in
    int *bitid,         //in
    bool *parent_check  //out
  ) {
  unsigned bitid_u  = static_cast<unsigned>(*bitid);
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    *parent_check = tree->block_is_parent(bitid_u);
  }
}

/** Wrapper function for MortonTree's identify */
extern "C" void bittree_identify(
    bool *updated,      //in
    int *lev,           //inout (0-based)
    int *ijk,           //inout
    int *mort,          //out
    int *bitid          //out
  ) {
  unsigned coord[BTDIM];
  for(unsigned d=0; d < BTDIM; d++)
    coord[d] = static_cast<unsigned>( ijk[d]);
  unsigned lev_u = static_cast<unsigned>(*lev);

  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    if(tree->inside(lev_u, coord)) {
      MortonTree::Block b = tree->identify(lev_u, coord);

      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < BTDIM; d++)
        ijk[d] = static_cast<int>(b.coord[d]);
      *mort = static_cast<int>(b.mort);
      *bitid = static_cast<int>(b.id);
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < BTDIM; d++)
        ijk[d] = -1;
      *mort = -1;
      *bitid = -1;
    }

  }
}

/** Wrapper function for MortonTree's locate */
extern "C" void bittree_locate(
    bool *updated,      //in
    int *bitid,         //in
    int *lev,           //out (0-based)
    int *ijk,           //out
    int *mort          //out
  ) {
  unsigned bitid_u  = static_cast<unsigned>(*bitid);

  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    if(bitid_u < tree->id_upper_bound() ) {
      MortonTree::Block b = tree->locate(bitid_u);

      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < BTDIM; d++)
        ijk[d] = static_cast<int>(b.coord[d]);
      *mort = static_cast<int>(b.mort);
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < BTDIM; d++)
        ijk[d] = -1;
      *mort = -1;
    }
  }
}

/** Get id0 */
extern "C" void bittree_get_id0(
    bool *updated,      //in
    int *idout          //out
  ) {
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    *idout = static_cast<int>(tree->level_id0(0));
  }
}

/** Wrapper function for TheTree's get_level_id_limits */
extern "C" void bittree_level_bitid_limits(
    bool *updated, //in
    int *lev,      //in
    int *ids       //out
  ) {
  unsigned lev_u  = static_cast<unsigned>(*lev);
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
    ids[0] = static_cast<int>(tree->level_id0(lev_u));
    ids[1] = static_cast<int>(tree->level_id1(lev_u));
  }
}

/** Wrapper function for TheTree's get_bitid_list, which 
  * itself wraps MortonTree's bitid_list */
extern "C" void bittree_get_bitid_list(
    bool *updated,      //in
    int *mort_min,      //in, 0-based
    int *mort_max,      //in, 0-based
    int *idout          //out
  ) {
  unsigned mort_min_u  = static_cast<unsigned>(*mort_min);
  unsigned mort_max_u  = static_cast<unsigned>(*mort_max);
  if(!!the_tree) {
    auto tree = the_tree->getTree(*updated);
#ifndef BITTREE_SAFE
    tree->bitid_list(mort_min_u, mort_max_u, idout);
#else
    int outlist[(mort_max_u - mort_min_u)];
    tree->bitid_list(mort_min_u, mort_max_u, outlist);
    for (unsigned i=0;i<(mort_max_u - mort_min_u);i++)
      idout[i] = outlist[i];
#endif
  }
}

/** Wrapper function for refine_init */
extern "C" void bittree_refine_init() {
  if(!!the_tree)
    the_tree->refine_init();
}

/** Wrapper funciton for refine_mark */
extern "C" void bittree_refine_mark(
    int *bitid,        // in
    bool *value        // in
  ) {
  unsigned bitid_u  = static_cast<unsigned>(*bitid);
  if(!!the_tree)
    the_tree->refine_mark(bitid_u, *value);
}

/** Wrapper function for refine_reduce */
extern "C" void bittree_refine_reduce(int *comm_) {
  if(!!the_tree) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_);
    the_tree->refine_reduce(comm);
  }
}

/** Wrapper function for refine_reduce_and */
extern "C" void bittree_refine_reduce_and(int *comm_) {
  if(!!the_tree) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_);
    the_tree->refine_reduce_and(comm);
  }
}

/** Wrapper function for refine_update */
extern "C" void bittree_refine_update() {
  if(!!the_tree)
    the_tree->refine_update();
}

/** Wrapper function for refine_apply */
extern "C" void bittree_refine_apply() {
  if(!!the_tree)
    the_tree->refine_apply();
}

/** Wrapper function the print_2d */
extern "C" void bittree_print(int *datatype)
{
  unsigned dtype_u = static_cast<unsigned>(*datatype);
  if(!!the_tree) {
    std::cout << the_tree->slice_to_string(dtype_u);
  }
}
