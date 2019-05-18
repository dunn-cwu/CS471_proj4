set -e
mkdir -p results/populations

echo "--------------"
./build/release/cs471-proj4.out ./params/hs-4-1k.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/hs-4-10k.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/hs-4-100k.ini
echo "--------------"

echo "Harmony search tests finished."
