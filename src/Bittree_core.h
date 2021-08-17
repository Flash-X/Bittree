#ifndef BITTREE_CORE_H__
#define BITTREE_CORE_H__

#include "Bittree_bitarray.h"
#include "Bittree_mortontree.h"
#include "Bittree_ref.h"
#include "mpi.h"
#include <iostream>

using namespace BitTree;

typedef unsigned W;

namespace { // private globals
  struct TheTree_ {
    virtual ~TheTree_() = default;
    virtual unsigned level_count(bool updated) = 0;
    virtual unsigned block_count(bool updated) = 0;
    virtual unsigned leaf_count(bool updated) = 0;
    virtual unsigned delta_count() = 0;
    virtual bool check_refine_bit(int bitid) = 0;
    virtual bool is_parent(bool updated, int bitid) = 0;
    virtual void identify(bool updated, int *lev, int *ijk, int *mort, int *bitid) = 0;
    virtual void locate(bool updated, int bitid, int *lev, int *ijk, int *mort) = 0;
    virtual void get_id0(bool updated, int *out) = 0;
    virtual void get_level_id_limits(bool updated, int lev, int *ids) = 0;
    virtual void get_bitid_list(bool updated, int mort_min, int mort_max, int *out) = 0;
    virtual void refine_init() = 0;
    virtual void refine_mark(int bitid, bool value) = 0;
    virtual void refine_reduce(MPI_Comm comm) = 0;
    virtual void refine_reduce_and(MPI_Comm comm) = 0;
    virtual void refine_update() = 0;
    virtual void refine_apply() = 0;
    virtual void print_2d(int datatype=0) = 0;
  };
  
  template<unsigned ndim>
  class TheTree: public TheTree_ {
    Ref<MortonTree<ndim,W> > tree;               //Actual Bittree
    Ref<MortonTree<ndim,W> > tree_updated;       //Updated Bittree, before refinement is applied
    Ref<BitArray<W> > refine_delta;              //(De)refinement flags for blocks
    bool is_reduced;                             //Flag to track whether refine_delta is up to date across processors
    bool is_updated;                             //Flag to track whether tree_updated matches latest refine_delta
    bool in_refine;                              //If in_refine=false, tree_updated and refine_delta should not exist 
  public:
    TheTree(const unsigned top[], const bool includes[]);
    unsigned level_count(bool updated);
    unsigned block_count(bool updated);
    unsigned leaf_count(bool updated);
    unsigned delta_count();
    bool check_refine_bit(int bitid);
    bool is_parent(bool updated, int bitid);
    void identify(bool updated, int *lev, int *ijk, int *mort, int *bitid);
    void locate(bool updated, int bitid, int *lev, int *ijk, int *mort);
    void get_id0(bool updated, int *out);
    void get_level_id_limits(bool updated, int lev, int *ids);
    void get_bitid_list(bool updated, int mort_min, int mort_max, int *out);
    void refine_init();
    void refine_mark(int bitid, bool value);
    void refine_reduce(MPI_Comm comm);
    void refine_reduce_and(MPI_Comm comm);
    void refine_update();
    void refine_apply();
    void print_2d(int datatype=0);
  };
  
  Ref<TheTree_> the_tree;
}

/** Checks if the_tree has been created */
extern "C" bool bittree_initialized();

/** Essentially a wrapper for TheTree's constructor */
extern "C" void bittree_init(
    int *ndim,      // in
    int topsize[],  // in
    bool includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
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

/** Wrapper function the print_2d */
extern "C" void bittree_print_2d(int *datatype=0);

#endif
