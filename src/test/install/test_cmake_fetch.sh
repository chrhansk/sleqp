set -e

cd cmake_fetch

mkdir -p build
cd build
cmake ..
cmake --build .
ctest
