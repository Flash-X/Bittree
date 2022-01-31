#ifndef BTUNIT_H__
#define BTUNIT_H__

#include <Bittree_Amr.h>
#include "constants.h"

using namespace bittree;

class btUnit {
  // Functions
  public:
    static void btErrorEst( std::shared_ptr<BittreeAmr> mesh );
    static void btRefineInitialize( std::shared_ptr<BittreeAmr> mesh );
    static void btRefineFinalize( std::shared_ptr<BittreeAmr> mesh );

  private:
    static void btCheckRefine( std::shared_ptr<BittreeAmr> mesh );
    static void btCheckDerefine( std::shared_ptr<BittreeAmr> mesh );
    static std::vector<int> calcNeighIntCoords(unsigned lev, unsigned* lcoord, int* gCell, std::shared_ptr<BittreeAmr> mesh);

  // cache of mesh data
  private:
    static std::vector<bool> refine;
    static std::vector<bool> derefine;
    static std::vector<std::vector<unsigned>> lcoord;
    static std::vector<unsigned> lev;
    static std::vector<bool> is_par;
    static std::vector<unsigned> bitid;
    static std::vector<std::vector<double>> error;
    static std::vector<std::vector<double>> error_par;

};


#endif
