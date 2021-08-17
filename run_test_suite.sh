#!/bin/bash

python setup.py test -d 1 --build build_1d --debug
python setup.py test -d 2 --build build_2d --debug
python setup.py test -d 3 --build build_3d --debug

cd build_1d
make
cd ..
cd build_2d
make
cd ..
cd build_3d
make
cd ..

cd build_1d
make test
exit_code_1=$?
cd ..
cd build_2d
make test
exit_code_2=$?
cd ..
cd build_3d
make test
exit_code_3=$?
cd ..

if [ $exit_code_1 -eq 0 ] && [ $exit_code_2 -eq 0 ] && [ $exit_code_3 -eq 0 ];
then
  echo SUCCESS!
  rm -r build_1d
  rm -r build_2d
  rm -r build_3d
else
  if [ $exit_code_1 -ne 0 ];
  then
    echo "Error: 1D test failed in execution"
  fi
  if [ $exit_code_2 -ne 0 ];
  then
    echo "Error: 2D test failed in execution"
  fi
  if [ $exit_code_3 -ne 0 ];
  then
    echo "Error: 3D test failed in execution"
  fi
  echo FAILURE!
fi
