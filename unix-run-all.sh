set -e
mkdir -p results/populations

./unix-run-pso.sh

echo "=============================================="

./unix-run-fa.sh

echo "=============================================="

./unix-run-hs.sh

echo "=============================================="
echo "All tests ran"
