set -e
mkdir -p results/populations

echo "--------------"
./build/release/cs471-proj4.out ./params/pso-1.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/pso-2.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/pso-3.ini
echo "--------------"
./build/release/cs471-proj4.out ./params/pso-4.ini
echo "--------------"

echo "Particle Swarm Tests Finished"
