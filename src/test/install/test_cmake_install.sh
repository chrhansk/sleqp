set -e

cd cmake_install/

mkdir -p build
cd build
cmake ..
cmake --build .
ctest
