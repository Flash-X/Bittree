#!/usr/bin/env python3

# Copyright 2022 UChicago Argonne, LLC and contributors
#
# Licensed under the Apache License, Version 2.0 (the "License"); 
# you may not use this file except in compliance with the License. 
#  
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, 
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and 
# limitations under the License.


# Users should invoke this script with a setup line similar to one of the following:
# `python setup.py test -d 2`
# `python setup.py library -d 2`
#
# To make the test, cd into the build directory and run `make` or `make all`. Then, the test can be run with
# `make test` and the code coverage report can be generated with `make coverage`.
#
# To get a summary of all the command line options, run `python setup.py --help`.

import argparse, sys, os, shutil
from subprocess import check_output

def main():
    # Parse command line args
    parser = argparse.ArgumentParser(description='Bittree Setup Tool')
    parser.add_argument('test',type=str,help='testName or library')
    parser.add_argument('--build','-b',type=str,default='build',help='Build directory.')
    parser.add_argument('--dim','-d',type=int,help='Dimensionality.')
    parser.add_argument('--debug',action="store_true",help='Set up in debug mode.')
    parser.add_argument('--coverage','-c',action="store_true",help='Enable code coverage.')
    parser.add_argument('--prefix',type=str,help='Where to install library.')
    args = parser.parse_args()

    print("Bittree setup")
    print("---------------------------")

    # Setup.py is located in the repository root directory.
    homeDir = os.path.dirname(os.path.abspath(sys.argv[0]))

    # Make build directory in root directory. Delete it first if it already exists.
    buildDir = os.path.join( homeDir, args.build)
    print("Creating build directory: "+args.build)
    if os.path.isdir(buildDir):
        shutil.rmtree(buildDir)
    os.makedirs(buildDir)

    # Link main makefile
    print("Linking Makefile")
    mainMakefile = os.path.join(homeDir,'Makefile')
    os.symlink(mainMakefile,os.path.join(buildDir,'Makefile'))

    # Link makefiles parts from site and src
    print("Linking Makefile.base")
    srcMakefile = os.path.join(homeDir,'src','Makefile.base')
    os.symlink(srcMakefile,os.path.join(buildDir,'Makefile.base'))

    siteMakefile = os.path.join(homeDir,'Makefile.site')
    if not os.path.isfile(siteMakefile):
        raise ValueError("Site Makefile not found in site directory")
    print("Linking Makefile.site")
    os.symlink(siteMakefile,os.path.join(buildDir,'Makefile.site'))

    # Find test directory
    testDir = os.path.join(homeDir,args.test)

    if (args.test == 'library'):
        if not args.prefix:
            raise ValueError("Need to supply prefix if building library!")
    else:
        if not os.path.isdir(testDir):
            raise ValueError("Test directory not found")

        # Get test makefile
        print("Linking Makefile.test from test: "+args.test)
        testMakefile = os.path.join(testDir,'Makefile.test')
        if not os.path.isfile(testMakefile):
            raise ValueError("Test Makefile not found in test dir")
        os.symlink(testMakefile,os.path.join(buildDir,'Makefile.test'))

    # Write Makefile.setup in builddir
    print("Writing Makefile.setup")
    setupMakefile = os.path.join(buildDir,'Makefile.setup')
    with open(setupMakefile,'w') as f:
        f.write("BASEDIR = {}\n".format(homeDir))
        f.write("BUILDDIR = $(BASEDIR)/{}\n".format(args.build))
        if args.debug:
            f.write("DEBUG = true\n")
        else:
            f.write("DEBUG = false\n")

        if args.coverage:
            f.write("CODECOVERAGE = true\n")
        else:
            f.write("CODECOVERAGE = false\n")

        f.write("BTDIM = {}\n".format(args.dim))

        f.write("\n")
        if args.test == 'library':
            f.write("LIBONLY = True\n")
            f.write("LIB_BITTREE_PREFIX = {}\n".format(args.prefix))
        else:
            f.write("# Leave blank if building a test\n")
            f.write("LIBONLY = \n")
            f.write("# Leave blank if not linking a prebuilt library!\n")
            f.write("LINKLIB = \n")
            f.write("# Should be current dir (i.e. `.`) if not linking prebuilt library\n")
            f.write("LIB_BITTREE = .\n")

    # Write Bittree_constants.h in builddir
    print("Writing Bittree_constants.h")
    constantsFile = os.path.join(buildDir,'Bittree_constants.h')
    with open(constantsFile,'w') as f:
        f.write("#ifndef BITTREE_CONSTANTS_H__\n#define BITTREE_CONSTANTS_H__\n\n")

        f.write("#define BTDIM       {}\n".format(args.dim))

        f.write("#endif\n")


    # Write the setup logfile
    print("Writing setup.log")
    logfileName = os.path.join(buildDir,"setup.log")
    with open(logfileName,'w') as f:
        f.write('Setup command line: \n')
        f.write(os.path.abspath(sys.argv[0]))
        f.write(' ')
        f.write(' '.join(sys.argv[1:]))
        f.write('\n\n\n')

        f.write('Build directory: \n')
        f.write(os.path.abspath(buildDir) )
        f.write('\n\n')

        f.write('Path to linked files:\n')
        f.write('Makefile --> {}\n'.format(os.path.abspath(mainMakefile)))
        f.write('Makefile.base --> {}\n'.format(os.path.abspath(srcMakefile)))
        f.write('Makefile.site --> {}\n'.format(os.path.abspath(siteMakefile)))
        if(args.test != 'library'):
            f.write('Makefile.test --> {}\n'.format(os.path.abspath(testMakefile)))
        f.write('\n')

        f.write('Contents of Makefile.setup:\n')
        f.write('----------------------------\n')
        with open(setupMakefile,'r') as fread:
            for line in fread:
                f.write(line)
        f.write('----------------------------\n\n')

        f.write('Contents of Bittree_constants.h:\n')
        f.write('----------------------------\n')
        with open(constantsFile,'r') as fread:
            for line in fread:
                f.write(line)
        f.write('----------------------------\n\n')

        f.write('Repository status:\n')
        f.write('----------------------------\n')
        f.write( check_output(['git','status']).decode('utf-8') )
        f.write('----------------------------\n\n')

        f.write('Most recent commits:\n')
        f.write('----------------------------\n')
        f.write( check_output(['git','log','--max-count','5']).decode('utf-8') )
        f.write('----------------------------\n')

    print("Successfully set up build directory!")


if __name__ == '__main__':
    main()
