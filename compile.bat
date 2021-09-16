:: Set the variables

set libs-path=C:\Local\LibsVs2017Py3
set visual-studio=C:\Program Files (x86)\Microsoft Visual Studio

:: Visual Studio shell

call "%visual-studio%\2017\Community\VC\Auxiliary\Build\vcvars64.bat"
call "%libs-path%\LibsVS2017Py3.cmd"

:: Creates the folders

rd /s /q build
mkdir build
cd build

:: Compile with CMake

cmake -A x64 -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release --target install

:: Make tests

ctest -j 8 -C Release -R PFEM3D_Metafor