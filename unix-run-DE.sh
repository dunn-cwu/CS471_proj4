set -e
mkdir -p results/diffevo
./build/release/cs471_proj3.out ./params/DE_Strat1.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat2.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat3.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat4.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat5.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat6.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat7.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat8.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat9.ini
echo "--------------"
./build/release/cs471_proj3.out ./params/DE_Strat10.ini
echo "--------------"
echo All tests ran

