set -e
mkdir -p results/populations

echo "--------------"
./build/release/cs471-proj4.out ./params/fa-1.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/fa-2.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/fa-3.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/fa-4.ini
echo "--------------"

echo "Firefly tests finished"
