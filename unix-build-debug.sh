set -e
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
make
mkdir -p debug
mv ./cs471_proj3.out debug/cs471_proj3.out
echo Program binary moved to build/debug
