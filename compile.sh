rm -rf build && mkdir build && cd build
cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=OFF ..
make -j8 install
