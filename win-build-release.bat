mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -G "MinGW Makefiles" ..
mingw32-make
mkdir release
move /y cs471_proj1.out.exe release\cs471_proj1.out.exe
cmd /k
