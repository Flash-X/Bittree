#include <Bittree_BittreeAmr.h>
#include <iostream>
#include <mpi.h>

using namespace bittree;


// Given a list of leaf blocks on a certain level, mark them for refinement.
// The points in the list are bottom left corner integer indices.
void mark_leaves(std::shared_ptr<BittreeAmr> mesh,
                   unsigned lev,
                   std::vector<std::vector<unsigned>> pts) {

    // Grab the pre-refinement tree
    auto tree0 = mesh->getTree(false);

    // For each point given, call MortontTree::identify to get the bitid,
    // then mark for nodetype change with BittreeAmr::refine_mark.
    for(unsigned i=0; i<pts.size(); ++i) {
        auto b = tree0->identify(lev, pts[i].data());
        if(lev!=b.level) continue; //the block does not exist
        mesh->refine_mark(b.id, true);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init ( &argc, &argv );
    MPI_Comm comm = MPI_COMM_WORLD;

    // Initialize mesh with 1 block on the top level.
    int top[BTDIM];
    for(unsigned d=0; d<BTDIM; ++d) top[d] = 1;
    int includes[1] = {1};
    auto mesh = std::make_shared<BittreeAmr>(top,includes);

    // First refinement step: create level 1 
    {
        // Call mesh->refine_init as the first step.
        mesh->refine_init();

        // This calls mesh->refine_mark on each block in the list.
        // For this first step, we're marking the only block that exists.
        std::vector<std::vector<unsigned>> pts = {{0,0}};
        mark_leaves(mesh,0,pts);

        // This reduces by ORing, so the list to be refined is the union of
        // blocks marked on each rank.
        mesh->refine_reduce(comm);

        // After all the tagging, a routine should be inserted here which checks
        // that all adjacent leaves differ by at most 1 level of refinement.

        // Generate the updated tree.
        mesh->refine_update();
        // At this point, the old and new tree coexist.

        std::cout << "Here's bittree during the first refinement step:\n";
        std::cout << "------------------------------------------------\n";
        std::cout << mesh->slice_to_string(0);

        // Destroy the old tree and exit refinement step.
        mesh->refine_apply();
    }

    // Create level 2
    {
        mesh->refine_init();
        std::vector<std::vector<unsigned>> pts = {{0,0},
                                                  {0,1},
                                                  {1,1}};
        mark_leaves(mesh,1,pts);
        mesh->refine_reduce(comm);
        mesh->refine_update();
        mesh->refine_apply();
    }

    // Create level 3
    {
        mesh->refine_init();
        std::vector<std::vector<unsigned>> pts = {{0,3},
                                                  {1,3}};
        mark_leaves(mesh,2,pts);
        mesh->refine_reduce(comm);
        mesh->refine_update();
        mesh->refine_apply();
    }

    // Create level 4
    {
        mesh->refine_init();
        std::vector<std::vector<unsigned>> pts = {{2,7}};
        mark_leaves(mesh,3,pts);
        mesh->refine_reduce(comm);
        mesh->refine_update();
        mesh->refine_apply();
    }

    std::cout << "Here's Bittree after all refinment steps:\n";
    std::cout << "-----------------------------------------\n";
    std::cout << mesh->slice_to_string(0);

    MPI_Finalize();
}
