#ifndef BTFUNC_H__
#define BTFUNC_H__

#include <Bittree_Amr.h>
#include "constants.h"

using namespace bittree;

class btUnit {
    public:


    // Functions
    static void btRefineInitialize( std::shared_ptr<BittreeAmr> mesh );
    static void btRefineFinalize( std::shared_ptr<BittreeAmr> mesh );


    private:
    // cache of mesh data
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
