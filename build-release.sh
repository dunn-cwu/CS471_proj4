cd build
cmake -DCMAKE_BUILD_TYPE=Release .
make
mkdir release
mv ./cs471_proj1.out release/cs471_proj1.out
