# Create the folders

rm -rf build
mkdir build
cd build

# Compile with CMake

cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release -DWITH_MPI=0 ..
make -j8 install

# Test the code

ctest -j4 -C Release -R PFEM3D_Metafor
