#include "Bittree_Amr.h"

using namespace bittree;

std::shared_ptr<MortonTree> BittreeAmr::getTree(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_;
  }
  else {
    return tree_;
  }
}


/** Constructor for BittreeAmr */
BittreeAmr::BittreeAmr(const unsigned top[], const bool includes[]):
  tree_(std::make_shared<MortonTree>(top, includes)),
  is_reduced_(false),
  is_updated_(false),
  in_refine_(false)  {
}

/** Check level count of Bittree. Wrapper for MortonTree's levels() */
unsigned BittreeAmr::level_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->levels();
  }
  else
    return tree_->levels();
}

/** Check block count of Bittree. Wrapper for MortonTree's blocks() */
unsigned BittreeAmr::block_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->blocks();
  }
  else
    return tree_->blocks();
}

/** Check block count of Bittree. Wrapper for MortonTree's leaves() */
unsigned BittreeAmr::leaf_count(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_->leaves();
  }
  else
    return tree_->leaves();
}
/** Check number of blocks marked for nodetype change */
unsigned BittreeAmr::delta_count() {
  if(in_refine_) return refine_delta_->count();
  else return 0;
}

/** Check refinement bit. Wrapper for BitArray's get() */
bool BittreeAmr::check_refine_bit(unsigned bitid) {
  if(in_refine_) return bool(refine_delta_->get(bitid));
  else return 0;
}


/** Wrapper for MortonTree's block_is_parent */ 
bool BittreeAmr::is_parent(bool updated, unsigned bitid) {
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
void BittreeAmr::identify(
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
void BittreeAmr::locate(
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

void BittreeAmr::get_id0(
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

void BittreeAmr::get_level_id_limits(
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
void BittreeAmr::get_bitid_list(
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
void BittreeAmr::refine_init() {
  unsigned nbits = tree_->id_upper_bound();
  refine_delta_ = std::make_shared<BitArray>(nbits);
  refine_delta_->fill(false);
  is_reduced_ = true;
  is_updated_ = false;
  in_refine_ = true;
}

/** Mark a bit on refine_delta_ */
void BittreeAmr::refine_mark(
    unsigned bitid,   // in
    bool value   // in
  ) {
  refine_delta_->set(bitid, value);
  is_reduced_ = false;
  is_updated_ = false;
}


/** Reduce refine_delta_ across all processors by ORing. This means
 *  any blocks marked on one processor will be marked on all. */
void BittreeAmr::refine_reduce(MPI_Comm comm) {
  int count = static_cast<int>(refine_delta_->word_count());
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta_->word_buf(),
    count,
    MPI_UNSIGNED,
    MPI_BOR,
    comm
  );
  is_reduced_ = true;
}

/** Reduce refine_delta_ across all processors by ANDing. This means
 *  any blocks unmarked on one processor will be unmarked on all. */
void BittreeAmr::refine_reduce_and(MPI_Comm comm) {
  int count = static_cast<int>(refine_delta_->word_count());
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta_->word_buf(),
    count,
    MPI_UNSIGNED,
    MPI_BAND,
    comm
  );
  is_reduced_ = true;
}

/** Generates the updated tree from the original tree + refine_delta_.
  * After call, both original and updated exist simultaneously. */
void BittreeAmr::refine_update() {
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
void BittreeAmr::refine_apply() {
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
void BittreeAmr::print_2d(unsigned datatype) {
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


