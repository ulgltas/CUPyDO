#!/bin/bash
# Configure script for CUPyDO
# External solver dir should place next to CUPyDO dir

. ~/.profile

export MTF_RUN="../Metafor"
export SU2_RUN="../SU2/bin"
# PATH
export PATH=$PATH:${SU2_RUN}
# LIB
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MTF_RUN}
# PYTHON
export PYTHONPATH="${PYTHONPATH}:${SU2_RUN}"
export PYTHONPATH="${PYTHONPATH}:../waves"
export PYTHONPATH="${PYTHONPATH}:../waves/build/bin"
export PYTHONPATH="${PYTHONPATH}:../NativeSolid/bin"
export PYTHONPATH=${PYTHONPATH}:${MTF_RUN}

echo ${PATH}
echo ${LD_LIBRARY_PATH}
echo ${PYTHON_PATH}

