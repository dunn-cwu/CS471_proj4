set -e
build/release/cs471_proj2.out params/param-blindsearch-10dim.ini
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-20dim.ini
echo "--------------"
build/release/cs471_proj2.out params/param-blindsearch-30dim.ini
echo "--------------"
echo All tests ran

