#include "Bittree_core.h"

using namespace bittree;


/** Constructor for TheTree */
TheTree::TheTree(const unsigned top[], const bool includes[]):
  tree_(std::make_shared<MortonTree>(top, includes)),
  is_reduced_(false),
  is_updated_(false),
  in_refine_(false)  {
}

/** Check level count of Bittree. Wrapper for MortonTree's levels() */
unsigned TheTree::level_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->levels();
  }
  else
    return tree_->levels();
}

/** Check block count of Bittree. Wrapper for MortonTree's blocks() */
unsigned TheTree::block_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->blocks();
  }
  else
    return tree_->blocks();
}

/** Check block count of Bittree. Wrapper for MortonTree's leaves() */
unsigned TheTree::leaf_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->leaves();
  }
  else
    return tree_->leaves();
}
/** Check number of blocks marked for nodetype change */
unsigned TheTree::delta_count() {
  if(in_refine_) return refine_delta_->count();
  else return 0;
}

/** Check refinement bit. Wrapper for BitArray's get() */
bool TheTree::check_refine_bit(unsigned bitid) {
  if(in_refine_) return bool(refine_delta_->get(bitid));
  else return 0;
}


/** Wrapper for MortonTree's block_is_parent */ 
bool TheTree::is_parent(bool updated, unsigned bitid) {
  if (updated && in_refine_){
    if (not is_updated_) refine_update();
    return tree_updated_->block_is_parent(bitid);
  }
  else
    return tree_->block_is_parent(bitid);
}

/** Identify a block based on its integer coordinates and level of refinement. 
 *  If no block exists on the level specified, Return the block on the finest level
 *  possible. */
void TheTree::identify(
    bool updated,        //in
    int *lev,            //inout (0-based)
    int *ijk,            //inout 
    int *mort,           //out
    int *bitid           //out
    ) {
  unsigned coord[NDIM];
  for(unsigned d=0; d < NDIM; d++)
    coord[d] = static_cast<unsigned>( ijk[d]);
  unsigned lev_u = static_cast<unsigned>(*lev);

  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    if(tree_updated_->inside(lev_u, coord)) {
      MortonTree::Block b = tree_updated_->identify(lev_u, coord);
      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = static_cast<int>(b.coord[d]);
      *mort = static_cast<int>(b.mort);
      *bitid = static_cast<int>(b.id);
    }
    else {
      *lev = -1;
      *mort = -1;
      *bitid = -1;
    }
  }
  else {
    if(tree_->inside(lev_u, coord)) {
      MortonTree::Block b = tree_->identify(lev_u, coord);
      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = static_cast<int>(b.coord[d]);
      *mort = static_cast<int>(b.mort);
      *bitid = static_cast<int>(b.id);
    }
    else {
      *lev = -1;
      *mort = -1;
      *bitid = -1;
    }
  }
}

/** Return location information about a block based on its bitid number. (Location in Bittree) 
 *   */
void TheTree::locate(
    bool updated,        //in
    unsigned bitid,           //in
    int *lev,            //out (0-based)
    int *ijk,            //out
    int *mort            //out
    ) {

  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    if(bitid < tree_updated_->id_upper_bound() ) {
      MortonTree::Block b = tree_updated_->locate(bitid);
      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = static_cast<int>(b.coord[d]);
      *mort = static_cast<int>(b.mort);
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = -1;
      *mort = -1;
    }
  }
  else {
    if(bitid < tree_->id_upper_bound() ) {
      MortonTree::Block b = tree_->locate(bitid);
      *lev = static_cast<int>(b.level);
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = static_cast<int>( b.coord[d]);
      *mort = static_cast<int>(b.mort);
    }
    else {
      *lev = -1;
      for(unsigned d=0; d < NDIM; d++)
        ijk[d] = -1;
      *mort = -1;
    }
  }
}

void TheTree::get_id0(
    bool updated,       //in
    int *out      //out
    ) {
  if(updated && in_refine_){
    if(not is_updated_) refine_update();
    out[0] = static_cast<int>(tree_updated_->level_id0(0));
  }
  else{
    out[0] = static_cast<int>(tree_->level_id0(0));
  }
}

void TheTree::get_level_id_limits(
    bool updated,   //in
    unsigned lev,        //in
    int *ids        //out
    ) {
  if(updated && in_refine_){
    if(not is_updated_) refine_update();
    tree_updated_->level_ids(lev, ids);
  }
  else{
    tree_->level_ids(lev, ids);
  }
}

/** Wrapper for MortonTree's bitid_list */
void TheTree::get_bitid_list(
    bool updated,       //in
    unsigned mort_min,  //in
    unsigned mort_max,  //in
    int *out      //out
    ) {
#ifndef BITTREE_SAFE
  if(updated && in_refine_){
    if(not is_updated_) refine_update();
    tree_updated_->bitid_list(mort_min, mort_max, out);
  }
  else{
    tree_->bitid_list(mort_min, mort_max, out);
  }
#else
  int outlist[(mort_max - mort_min)];
  if(updated && in_refine_){
    if(not is_updated_) refine_update();
    tree_updated_->bitid_list(mort_min, mort_max, outlist);
  }
  else{
    tree_->bitid_list(mort_min, mort_max, outlist);
  }
  for (unsigned i=0;i<(mort_max - mort_min);i++)
    out[i] = outlist[i];
#endif
}


/** Creates refine_delta_, and initializes all values to False.
  * First step of refinement. */
void TheTree::refine_init() {
  unsigned nbits = tree_->id_upper_bound();
  refine_delta_ = std::make_shared<BitArray>(nbits);
  refine_delta_->fill(false);
  is_reduced_ = true;
  is_updated_ = false;
  in_refine_ = true;
}

/** Mark a bit on refine_delta_ */
void TheTree::refine_mark(
    unsigned bitid,   // in
    bool value   // in
  ) {
  refine_delta_->set(bitid, value);
  is_reduced_ = false;
  is_updated_ = false;
}


/** Reduce refine_delta_ across all processors by ORing. This means
 *  any blocks marked on one processor will be marked on all. */
void TheTree::refine_reduce(MPI_Comm comm) {
  int count = static_cast<int>(refine_delta_->word_count());
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta_->word_buf().data(),
    count,
    MPI_UNSIGNED,
    MPI_BOR,
    comm
  );
  is_reduced_ = true;
}

/** Reduce refine_delta_ across all processors by ANDing. This means
 *  any blocks unmarked on one processor will be unmarked on all. */
void TheTree::refine_reduce_and(MPI_Comm comm) {
  int count = static_cast<int>(refine_delta_->word_count());
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta_->word_buf().data(),
    count,
    MPI_UNSIGNED,
    MPI_BAND,
    comm
  );
  is_reduced_ = true;
}

/** Generates the updated tree from the original tree + refine_delta_.
  * After call, both original and updated exist simultaneously. */
void TheTree::refine_update() {
  if (not is_reduced_) {
    std::cout << "Bittree updating before reducing. Possible error." << std::endl;
  }
  //std::cout << "refine_delta_->length() = " << refine_delta_->length() << std::endl;
  //for (unsigned j=0; j<refine_delta_->length() ; j++){
  //  std::cout << j << ": " << refine_delta_->get(j) << ";  " ;
  //}
  //std::cout << std::endl;
  tree_updated_ = tree_->refine(refine_delta_);
  is_updated_ = true;
}

/** Makes the updated tree the original tree. Final step of refinement. */
void TheTree::refine_apply() {
  if (not is_updated_) {
    refine_update();
  }
  tree_ = tree_updated_;
  refine_delta_ = nullptr;
  tree_updated_ = nullptr;
  in_refine_ = false;
}

/** Wrapper function to MortonTree's print_if_2d, which print a nice 
  * representation of the (2d) Bittree and refine_delta_. 
  * If tree has been updated, print both original and updated version. 
  * Can be passed a datatype to change what number prints at each block loc. 
  * (0=bitid, 1=morton number, 2=parentage) */
void TheTree::print_2d(unsigned datatype) {
  if (NDIM==2) { 
    switch (datatype) {
      case 0: 
        std::cout << "printing original tree (datatype=bitid): " <<std::endl;
        break;
      case 1:
        std::cout << "printing original tree (datatype=mort): " <<std::endl;
        break;
      case 2:
        std::cout << "printing original tree (datatype=parent): " <<std::endl;
        break;
    }
    tree_->print_if_2d(datatype);
    if(in_refine_) {
    std::cout << "printing refine_delta_ (indexed by bitid):" <<std::endl;
    for (unsigned j=0; j<refine_delta_->length() ; j++){
      std::cout << j << ": " << refine_delta_->get(j) << ";  " ;
    }
    }
    std::cout << std::endl;
    if (in_refine_ && is_updated_) {
      switch (datatype) {
      case 0: 
        std::cout << "printing updated tree (datatype=bitid): " <<std::endl;
        break;
      case 1:
        std::cout << "printing updated tree (datatype=mort): " <<std::endl;
        break;
      case 2:
        std::cout << "printing updated tree (datatype=parent): " <<std::endl;
        break;
      }
      tree_updated_->print_if_2d(datatype);
    }
  }
  else {
    std::cout << "Error: tried to print 2d Bittree but NDIM is not 2!" <<std::endl;
  }
}


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
 
  the_tree = std::make_shared<TheTree>(top,includes);
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

