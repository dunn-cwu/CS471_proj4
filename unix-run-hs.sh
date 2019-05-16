set -e
mkdir -p results/populations

echo "--------------"
./build/release/cs471-proj4.out ./params/hs-1.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/hs-2.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/hs-3.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/hs-4.ini
echo "--------------"

echo "Harmony search tests finished."
