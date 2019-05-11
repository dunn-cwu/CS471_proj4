set -e
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
mkdir -p release
mv ./cs471_proj4.out release/cs471_proj4.out
echo Program binary moved to build/release
