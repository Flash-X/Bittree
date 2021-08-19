#include "Bittree_fi.h"

/** Checks if the_tree has been created */
extern "C" bool bittree_initialized() {
  return !!the_tree;
}

/** Essentially a wrapper for TheTree's constructor */
extern "C" void bittree_init(
    int topsize[],  // in
    bool includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
  ) {
  unsigned top[3];
  for(int d=0; d < NDIM; d++)
    top[d] = static_cast<unsigned>(topsize[d]);
 
  the_tree = std::make_shared<BittreeAmr>(top,includes);
}

/** Wrapper function for block_count */
extern "C" void bittree_level_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree)
    *count = static_cast<int>(the_tree->level_count(*updated));
}


/** Wrapper function for block_count */
extern "C" void bittree_block_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree)
    *count = static_cast<int>(the_tree->block_count(*updated));
}


/** Wrapper function for leaf_count */
extern "C" void bittree_leaf_count(
    bool *updated,     //in: boolean
    int *count         //out
  ) {
  if(!!the_tree)
    *count = static_cast<int>(the_tree->leaf_count(*updated));
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
  if(!!the_tree)
    *parent_check = the_tree->is_parent(*updated,bitid_u);
}

/** Wrapper function for TheTree's identify, which 
  * itself wraps MortonTree's identify */
extern "C" void bittree_identify(
    bool *updated,      //in
    int *lev,           //inout (0-based)
    int *ijk,           //inout
    int *mort,          //out
    int *bitid          //out
  ) {
  if(!!the_tree)
    the_tree->identify(*updated, lev, ijk, mort, bitid);
}

/** Wrapper function for TheTree's locate, which 
  * itself wraps MortonTree's locate */
extern "C" void bittree_locate(
    bool *updated,      //in
    int *bitid,         //in
    int *lev,           //out (0-based)
    int *ijk,           //out
    int *mort          //out
  ) {
  unsigned bitid_u  = static_cast<unsigned>(*bitid);
  if(!!the_tree)
    the_tree->locate(*updated, bitid_u, lev, ijk, mort);
}

/** Wrapper function for TheTree's get_id0 */
extern "C" void bittree_get_id0(
    bool *updated,      //in
    int *idout          //out
  ) {
  if(!!the_tree)
    the_tree->get_id0(*updated, idout);
}

/** Wrapper function for TheTree's get_level_id_limits */
extern "C" void bittree_level_bitid_limits(
    bool *updated, //in
    int *lev,      //in
    int *ids       //out
  ) {
  unsigned lev_u  = static_cast<unsigned>(*lev);
  if(!!the_tree)
    the_tree->get_level_id_limits(*updated, lev_u, ids);
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
  if(!!the_tree)
    the_tree->get_bitid_list(*updated, mort_min_u, mort_max_u, idout);
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
extern "C" void bittree_print_2d(int *datatype)
{
  unsigned dtype_u = static_cast<unsigned>(*datatype);
  if(!!the_tree)
    the_tree->print_2d(dtype_u) ;
}
