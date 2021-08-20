#include "btUnit.h"
#include <iostream>

// static variable initialization
std::vector<bool> btUnit::refine;
std::vector<bool> btUnit::derefine;
std::vector<std::vector<unsigned>> btUnit::lcoord;
std::vector<unsigned> btUnit::lev;
std::vector<unsigned> btUnit::bitid;
std::vector<bool> btUnit::is_par;
std::vector<std::vector<double>> btUnit::error;
std::vector<std::vector<double>> btUnit::error_par;


void btUnit::btRefineInitialize( std::shared_ptr<BittreeAmr> mesh ) {
    // Tree before refinement. With only one rank, lnblocks = nblocks.
    auto tree0 = mesh->getTree();
    unsigned lnblocks = tree0->blocks();

    // Initialize the mesh data caches to appropriate values
    refine    = std::vector<bool>(lnblocks, false);
    derefine  = std::vector<bool>(lnblocks, false);
    lcoord    = std::vector<std::vector<unsigned>>(
                      lnblocks, std::vector<unsigned>(NDIM,0u) );
    lev       = std::vector<unsigned>(lnblocks, 0 );
    is_par    = std::vector<bool>(lnblocks, false);
    bitid     = std::vector<unsigned>(lnblocks, 0 );
    error     = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );
    error_par = std::vector<std::vector<double>>(
                      lnblocks, std::vector<double>(NVARS,0.0) );


    // Loop over local blocks to cache metadata in order of morton curve.
    // In this case I used bittree functionality to actually get the
    // coordinates and level of the list of blocks, but
    // in AMReX this data can be gotten elsewhere obviously.
    unsigned id0 = tree0->level_id0(0);
    unsigned id1 = tree0->id_upper_bound();
    for( unsigned lb = id0; lb < id1; ++lb) {
        MortonTree::Block actual = tree0->locate(lb);

        // Bittree calcualtes blkID - the local ID the block would have
        // in paramesh
        MortonTree::Block b = tree0->identify( actual.level, actual.coord);
        unsigned blkID = b.mort;

        // Should probably check to make sure actual.lev == b.lev

        // Cache lev, lcoord, and is_par for later reference
        for(unsigned d=0; d<NDIM; ++d)
            lcoord[blkID][d] = actual.coord[d];
        lev[blkID] = actual.level;
        is_par[blkID] = tree0->block_is_parent(b.id);
        bitid[blkID] = b.id;

        // Estimate error
        unsigned error_calc_result;
        for(unsigned v=0; v<NVARS; ++v) {
            error_calc_result = 0; // replace with actual calculation
            error[blkID][v] = error_calc_result;
        }
    }


    // Exchange error between parents and children. On one rank, this is
    // just having each parent give its error to all of its children.
    for( unsigned v=0; v<NVARS; ++v) {
      for( unsigned lb = 0; lb<lnblocks; ++lb) { //lb now is local blk num
        if(is_par[lb]) {
          unsigned ch[3];
          unsigned lev_ch;
          unsigned coord_ch[NDIM];
          for(ch[2]=0; ch[2]<=K3D; ++ch[2]) {
          for(ch[1]=0; ch[1]<=K2D; ++ch[1]) {
          for(ch[0]=0; ch[0]<=K1D; ++ch[0]) {
              lev_ch = lev[lb] + 1;
              for(unsigned d=0; d<NDIM; ++d) {
                  coord_ch[d] = lcoord[lb][d]*2 + ch[d];
              }
              auto b = tree0->identify(lev_ch, coord_ch);

              error_par[b.mort][v] = error[lb][v];
          }}}
        }
      }
    }

    // Actual marking for refinement/derefinement
    unsigned refineCutoff, derefineCutoff;
    bool allVarsDeref;
    unsigned maxLevel = 4;
    unsigned minLevel = 1;

    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      allVarsDeref = true;
      for(unsigned v=0; v<NVARS; ++v) { 

        // These are application parameters, could be different per var
        refineCutoff = 1.0;
        derefineCutoff = 1.0; 

        // If block's error is too large for any single refinement variable,
        // then the block should be refined. The block's error is too small
        // for ALL of the variables, the block should be derefined.
        if( lev[lb]<maxLevel &&
            !is_par[lb] &&
            (error[lb][v]>refineCutoff || lev[lb]<minLevel)) {
          refine[lb] = true;
          derefine[lb] = false; 
        }

        if(refine[lb]) break; // no need to do other vars

        if( lev[lb] > minLevel  && !is_par[lb] ) {
          if((error[lb][v] <= derefineCutoff  &&
              error_par[lb][v] <= derefineCutoff && allVarsDeref)
             || lev[lb] > maxLevel ) {
            derefine[lb] = true;
          }
        }

        if( error[lb][v] > derefineCutoff ||
            error_par[lb][v] > derefineCutoff) {
           allVarsDeref = false;
           derefine[lb] = false;
        }
      }
    }

//---------------------------------------------------------------------
//--Initialize bittree refinement and mark leaves to be refined
    mesh->refine_init();

    for( unsigned lb = 0; lb<lnblocks; ++lb) {
      if (refine[lb])
        mesh->refine_mark(bitid[lb], true);
    }

    // check_refine

    // mark for derefinement

    // check_derefine

//---------------------------------------------------------------------
//--Generate updated Bittree
    mesh->refine_update();

    // sort (distribute over procs)

}


void btUnit::btRefineFinalize( std::shared_ptr<BittreeAmr> mesh ) {
    mesh->refine_apply();
    // verify tree?
}
