# Tutorial: How to prepare CMakeLists.txt for your application program

This directory contains an example of CMakeLists.txt for user application program using Rokko.

## How to compile and test

```
cmake -DROKKO_ROOT_DIR=/place/where/rokko/is/installed .
make
ctest
```

or

```
export ROKKO_ROOT=/place/where/rokko/is/installed
cmake .
make
ctest
```
