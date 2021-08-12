# What is Bittree?
From John Bachan's paper: "The bittree is a memory efficient data structure for storing the entire octree structure of the mesh on each node. It also yields an efficient algorithm for mapping from a block's coordinates in the domain to its index along the space filling Morton curve which is used to distribute blocks among processors. Together, these features enable determination of the processor to which a given block belongs without requiring any off-node communication."

# Bittree Unit Test
This repository contains the Bittree source code (in `src`) as well as a unit test built with GoogleTest. Before building the test, customize Makefile.site by filling in appropriate information for your system. Then in the repository root directory, run `make` to build the test. The test can be executed with `make test`. As long as `CODECOVERAGE=true` is specified in Makefile.site, a coverage report can be generated with `make coverage`.
