This is a tutorial for creating and refining a mesh. I arbitrarily defined an AMR mesh focused around a feature toward the top left. In the file `diagram.pdf` you can see the grid and the breakdown of each level. In `main.cpp` I just create a BittreeAmr object and refine it consecutively so it matches with the diagram, which you can confirm with the standard output.

Before compiling
================
1. Compile and install Bittree library as per top-level README.
2. Edit the first section of the Makefile here with appropriate paths.


After that you should be good to compile and run the program.
