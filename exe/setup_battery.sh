# Other code versions

export SU2=fix_wrap_strong
export PYBEAM=master
export GEOGEN=v2.0.0
export MODALI=v2.0
export VLM=v2.0
export PFEM=v1.27
export NATIVESOLID=v1.1
export DARTFLO=master

# Development code versions

export OOMETA=master
export OONDA=master
export PARASOLID=master
export LINUXBIN=master

export PFEM3D=v2.4.0
export CUPYDO=master

# Gitlab and GitHub access token

export TOKEN_LAB=
export TOKEN_HUB=
export LAB=username:${TOKEN_LAB}@
export HUB=username:${TOKEN_HUB}@

# Creates output folder

export OUTPUT=Battery
rm -rf ${OUTPUT} && mkdir ${OUTPUT}
export OUTPUT=$PWD/${OUTPUT}
cd ${OUTPUT}

# Shortcut for cloning a repository

Clone(){
    
    git clone -c advice.detachedHead=false \
    --recursive --branch $1 https://$2.git
}

# Download Modali on GitHub

Modali(){

    cd ${OUTPUT}
    Clone ${MODALI} github.com/ulgltas/modali
}

# Download NativeSolid on GitHub

NativeSolid(){

    cd ${OUTPUT}
    Clone ${NATIVESOLID} github.com/ulgltas/NativeSolid
    
    cd NativeSolid && mkdir build && cd build
    cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release ..
    make -j8
}

# Download PyBeam on GitHub

pyBeam(){
    
    cd ${OUTPUT}
    Clone ${PYBEAM} github.com/pyBeam/pyBeam

    cd pyBeam && rm -rf bin && mkdir build
    sed -i "s/dependency('python3')/dependency('python3-embed')/g" meson.build
    meson build --prefix=${OUTPUT}/pyBeam/build
    ninja -C build install
}

# Download VLM on GitHub

VLM(){

    cd ${OUTPUT}
    Clone ${VLM} github.com/ulgltas/VLM
    Clone ${GEOGEN} github.com/acrovato/geoGen
    
    cd VLM && mkdir build && cd build
    cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release -DPYTHON_WRAPPER=1 \
    -DCMAKE_INSTALL_PREFIX=bin ..
    make -j8 install
}

# Download SU2 on GitHub

SU2(){

    cd ${OUTPUT}
    Clone ${SU2} github.com/ulgltas/SU2
    cd SU2 && mkdir build
    
    python3 meson.py setup build \
    -Denable-pywrapper=true \
    -Denable-cgns=false \
    -Denable-tecio=false \
    -Dwith-mpi=disabled \
    -Denable-tests=false \
    -Denable-autodiff=true \
    --prefix=${OUTPUT}/SU2/build
    ./ninja -C build install
}

# Download PFEM3D on GitHub

PFEM3D(){

    cd ${OUTPUT}
    Clone ${PFEM3D} ${HUB}github.com/ImperatorS79/PFEM3D

    cd PFEM3D && mkdir build && cd build
    cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=1 -DUSE_SWIG=1 -DUSE_MKL=1 ..
    make -j8
}

# Download DARTFlo on GitLab

DARTFlo(){

    cd ${OUTPUT}
    Clone ${DARTFLO} ${LAB}gitlab.uliege.be/am-dept/dartflo
    
    cd dartflo && mkdir build && cd build
    cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=bin ..
    make -j8 install
}

# Download PFEM on GitLab

PFEM(){

    cd ${OUTPUT}
    Clone ${PFEM} ${LAB}gitlab.uliege.be/am-dept/PFEM

    cd PFEM && mkdir build && cd build
    cmake -Wno-dev -DPFEM_USE_MUMPS=1 \
    -DCMAKE_BUILD_TYPE=Release ..
    make -j8
}

# Download Metafor on GitLab

Metafor(){

    cd ${OUTPUT}
    mkdir Metafor && cd Metafor

    Clone ${OONDA} ${LAB}gitlab.uliege.be/am-dept/MN2L/oo_nda
    Clone ${OOMETA} ${LAB}gitlab.uliege.be/am-dept/MN2L/oo_meta
    Clone ${PARASOLID} ${LAB}gitlab.uliege.be/am-dept/MN2L/parasolid
    Clone ${LINUXBIN} ${LAB}gitlab.uliege.be/am-dept/linuxbin

    mkdir oo_metaB && cd oo_metaB
    cmake -Wno-dev -C ../oo_meta/CMake/vaillant.cmake \
    -DCMAKE_BUILD_TYPE=Release ../oo_meta
    make -j8
}

# Download CUPyDO on GitHub

CUPyDO(){

    cd ${OUTPUT}
    Clone ${CUPYDO} github.com/ulgltas/CUPyDO

    cd CUPyDO && mkdir build && cd build
    cmake -Wno-dev -DCMAKE_BUILD_TYPE=Release \
    -DWITH_MPI=0 -DCMAKE_INSTALL_PREFIX=bin ..
    make -j8 install
}

# Call the Functions

Modali && VLM && SU2 && DARTFlo
PFEM && pyBeam && NativeSolid
PFEM3D && Metafor && CUPyDO
