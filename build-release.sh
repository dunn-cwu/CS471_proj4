set -e
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
mkdir -p release
mv ./cs471_proj2.out release/cs471_proj2.out
echo Program binary moved to build/release
