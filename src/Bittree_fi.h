/*
   Licensed under the Apache License, Version 2.0 (the "License"); 
   you may not use this file except in compliance with the License. 
    
 
   Unless required by applicable law or agreed to in writing, software 
   distributed under the License is distributed on an "AS IS" BASIS, 
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
   See the License for the specific language governing permissions and 
   limitations under the License.
*/
#ifndef BITTREE_FI_H__
#define BITTREE_FI_H__

#include "Bittree_BittreeAmr.h"

using namespace bittree;

namespace {
    std::shared_ptr<BittreeAmr> the_tree;
}

/** Checks if the_tree has been created */
extern "C" bool bittree_initialized();

/** Essentially a wrapper for TheTree's constructor */
extern "C" void bittree_init(
    int topsize[],  // in
    int includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
  );

/** Wrapper function for block_count */
extern "C" void bittree_level_count(
    bool *updated,     //in: boolean
    int *count         //out
  );


/** Wrapper function for block_count */
extern "C" void bittree_block_count(
    bool *updated,     //in: boolean
    int *count         //out
  );


/** Wrapper function for leaf_count */
extern "C" void bittree_leaf_count(
    bool *updated,     //in: boolean
    int *count         //out
  );

/** Wrapper function for delta_count */
extern "C" void bittree_delta_count(
    int *count         //out
  );


/** Wrapper function for check_refine_bit */
extern "C" void bittree_check_refine_bit(
    const int *bitid,   //in
    bool *bit_check     //out
  );

/** Wrapper function for is_parent */
extern "C" void bittree_is_parent(
    bool *updated,      //in
    int *bitid,         //in
    bool *parent_check  //out
  );

/** Wrapper function for TheTree's identify, which 
  * itself wraps MortonTree's identify */
extern "C" void bittree_identify(
    bool *updated,      //in
    int *lev,           //inout (0-based)
    int *ijk,           //inout
    int *mort,          //out
    int *bitid          //out
  );

/** Wrapper function for TheTree's locate, which 
  * itself wraps MortonTree's locate */
extern "C" void bittree_locate(
    bool *updated,      //in
    int *bitid,         //in
    int *lev,           //out (0-based)
    int *ijk,           //out
    int *mort          //out
  );

/** Wrapper function for TheTree's get_id0 */
extern "C" void bittree_get_id0(
    bool *updated,      //in
    int *idout          //out
  );

/** Wrapper function for TheTree's get_level_id_limits */
extern "C" void bittree_level_bitid_limits(
    bool *updated, //in
    int *lev,      //in
    int *ids       //out
  );

/** Wrapper function for TheTree's get_bitid_list, which 
  * itself wraps MortonTree's bitid_list */
extern "C" void bittree_get_bitid_list(
    bool *updated,      //in
    int *mort_min,      //in, 0-based
    int *mort_max,      //in, 0-based
    int *idout          //out
  );

/** Wrapper function for refine_init */
extern "C" void bittree_refine_init();

/** Wrapper funciton for refine_mark */
extern "C" void bittree_refine_mark(
    int *bitid,        // in
    bool *value        // in
  );

/** Wrapper function for refine_reduce */
extern "C" void bittree_refine_reduce(int *comm_);

/** Wrapper function for refine_reduce_and */
extern "C" void bittree_refine_reduce_and(int *comm_);

/** Wrapper function for refine_update */
extern "C" void bittree_refine_update();

/** Wrapper function for refine_apply */
extern "C" void bittree_refine_apply();

/** print (slice=0) */
extern "C" void bittree_print(int *datatype=0);

#endif
