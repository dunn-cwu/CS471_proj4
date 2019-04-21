set -e
build/release/cs471_proj2.out params/param-localsearch-10dim.ini 0
echo "--------------"
build/release/cs471_proj2.out params/param-localsearch-20dim.ini 0
echo "--------------"
build/release/cs471_proj2.out params/param-localsearch-30dim.ini 0
echo "--------------"
echo All tests ran

