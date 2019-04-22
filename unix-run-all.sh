set -e
./build/release/cs471_proj2.out ./params/param-blindsearch-10dim.ini
echo "--------------"
./build/release/cs471_proj2.out ./params/param-blindsearch-20dim.ini
echo "--------------"
./build/release/cs471_proj2.out ./params/param-blindsearch-30dim.ini
echo "--------------"
./build/release/cs471_proj2.out ./params/param-localsearch-10dim.ini
echo "--------------"
./build/release/cs471_proj2.out ./params/param-localsearch-20dim.ini
echo "--------------"
./build/release/cs471_proj2.out ./params/param-localsearch-30dim.ini
echo "--------------"
echo All tests ran

