#include <iostream>
#include "bittree_core.cxx"

int main() {
    int ndim = 2;
    int top[2] = {2,3};
    bool includes[6];
    for(int i=0; i<6; i++) includes[i] = true;


    bool updated = false;
    int count = 0;
    bittree_init(&ndim, top, includes);

    bittree_block_count(&updated, &count);

    std::cout << "Count = " << count << std::endl;
}
