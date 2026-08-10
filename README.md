# What is Bittree?
From John Bachan's paper: "The bittree is a memory efficient data structure for storing the entire octree structure of the mesh on each node. It also yields an efficient algorithm for mapping from a block's coordinates in the domain to its index along the space filling Morton curve which is used to distribute blocks among processors. Together, these features enable determination of the processor to which a given block belongs without requiring any off-node communication."

# Bittree Unit Test
This repository contains the Bittree source code (in `src`) as well as a unit test built with GoogleTest. Before building the test, customize Makefile.site by filling in appropriate information for your system. Then in the repository root directory, run the run\_test\_suite.sh bash script which builds and executes 1D, 2D, and 3D versions of the test.

# Building Bittree Library
After customizing Makefile.site, use the following commands to build and install the Bittree library.

```
python setup.py library --dim N --prefix NameOfPrefix
cd build
make
make install
```

# Building with CMake

CMake can be used instead of `setup.py`, and needs no `Makefile.site`:

```
mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX="<prefix>" -DBITTREE_NDIM="2;3"
make
ctest
make install
```

`BITTREE_NDIM` is the list of dimensions to build (default `1;2;3`), and each dimension N is
installed under `<prefix>/Nd`. Other options are `BITTREE_BUILD_TESTS`,
`BITTREE_ENABLE_COVERAGE`, `BITTREE_WERROR`, and `BITTREE_ALT_MORTON_ORDER`.

Another CMake project uses an installed Bittree by pointing at the dimension it needs,

```
cmake -S . -B build -DBittree_ROOT=<prefix>/2d
```

and linking the imported target.

```
find_package(Bittree REQUIRED COMPONENTS 2d)
target_link_libraries(myapp PRIVATE Bittree::bittree)
```

`COMPONENTS 2d` is optional but recommended: the build then fails with a clear message if it is
pointed at a Bittree of the wrong dimension, instead of silently linking it.

# Bittree Tutorial

The Bittree examples in the `tutorial` directory requires the 2D library to be built first. Then go the `Makefile` and appropriately fill in the the top section. The test can be made with `make` and run with `make test`.
