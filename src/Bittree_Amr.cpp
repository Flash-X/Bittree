#include "Bittree_Amr.h"
#include <sstream>

namespace bittree {

/** Constructor for BittreeAmr */
BittreeAmr::BittreeAmr(const unsigned top[], const bool includes[]):
  tree_(std::make_shared<MortonTree>(top, includes)),
  is_reduced_(false),
  is_updated_(false),
  in_refine_(false)  {
}

/** Get shared_ptr to the actual Bittree.
  * If in the middle of refinement, can also get the updated tree.
  */
std::shared_ptr<MortonTree> BittreeAmr::getTree(bool updated) {
  if(updated && in_refine_) {
    if (not is_updated_) refine_update();
    return tree_updated_;
  }
  else {
    return tree_;
  }
}

/** Check number of blocks marked for nodetype change */
unsigned BittreeAmr::delta_count() const {
  if(in_refine_) return refine_delta_->count();
  else return 0;
}

/** Check refinement bit. Wrapper for BitArray's get() */
bool BittreeAmr::check_refine_bit(unsigned bitid) const {
  if(in_refine_) return bool(refine_delta_->get(bitid));
  else return 0;
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
  if(in_refine_) {
      refine_delta_->set(bitid, value);
  }
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

/** Wrapper function to MortonTree's print_slice, which print a nice
  * representation of the Bittree and refine_delta_.
  * If tree has been updated, print both original and updated version.
  * Can be passed a datatype to change what number prints at each block loc.
  * (0=bitid, 1=morton number, 2=parentage) */
std::string BittreeAmr::slice_to_string(unsigned datatype, unsigned slice) const {
  std::ostringstream buffer;
  if (datatype==0 || datatype==1 || datatype==2) {
    buffer << "printing original tree:\n";
    buffer << tree_->print_slice(datatype,slice) << "\n";

    if(in_refine_) {
    buffer << "printing refine_delta_ (indexed by bitid):\n";
    unsigned id0 = tree_->level_id0(0);
    for (unsigned j=id0; j<refine_delta_->length() ; j++){
      buffer << j << ": " << refine_delta_->get(j) << ";  " ;
    }
    buffer << "\n\n";
    }

    if (in_refine_ && is_updated_) {
      buffer << "printing updated tree:\n";
      buffer << tree_updated_->print_slice(datatype,slice) << "\n";
    }
  }
  else {
    buffer << "Error: datatype must be 0(bitid), 1(mort), or 2(parent)!\n";
  }
  return buffer.str();
}

}



