set -e
build/release/cs471_proj2.out params/param-blindsearch-10dim.ini 0
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-20dim.ini 0
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-30dim.ini 0
echo "--------------"
echo All tests ran

