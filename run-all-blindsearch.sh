set -e
build/release/cs471_proj2.out params/param-blindsearch-10dim.ini 2
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-20dim.ini 2
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-30dim.ini 2
echo "--------------"
echo All tests ran

