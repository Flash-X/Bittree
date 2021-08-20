#include <Bittree_Amr.h>
#include "btUnit.h"
#include "constants.h"

#include <iostream>

using namespace bittree;

// local global access to the mesh
namespace {
    std::shared_ptr<BittreeAmr> mesh;
}

int main() {
    unsigned top[NDIM];
    for(unsigned d=0; d<NDIM; ++d) top[d] = 1;
    bool includes[1] = {true};
    mesh = std::make_shared<BittreeAmr>(top,includes);

    btUnit::btRefineInitialize( mesh );
    std::cout << mesh->slice_to_string(0);
    btUnit::btRefineFinalize( mesh );

}
