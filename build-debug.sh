set -e
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../
make
mkdir -p debug
mv ./cs471_proj1.out debug/cs471_proj1.out
echo Program binary moved to build/debug
