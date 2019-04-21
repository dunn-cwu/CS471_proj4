set -e
build/release/cs471_proj2.out params/param-localsearch-10dim.ini 1
echo "--------------"
build/release/cs471_proj2.out params/param-localsearch-20dim.ini 1
echo "--------------"
build/release/cs471_proj2.out params/param-localsearch-30dim.ini 1
echo "--------------"
echo All tests ran

