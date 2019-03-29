cd build
cmake -DCMAKE_BUILD_TYPE=Debug .
make
mkdir debug
mv ./cs471_proj1.out debug/cs471_proj1.out
