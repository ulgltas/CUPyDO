#!/bin/bash
# Configure script for CUPyDO
# External solver dir should place next to CUPyDO dir

. ~/.profile

export TOP_DIR=$(cd "$(dirname "$0")/.."; pwd)
export MTF_RUN="${TOP_DIR}/Metafor"
export SU2_RUN="${TOP_DIR}/SU2/bin"
# PATH
export PATH=$PATH:${SU2_RUN}
# LIB
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MTF_RUN}
# PYTHON
export PYTHONPATH=${PYTHONPATH}:${SU2_RUN}
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/waves"
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/waves/build/bin"
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/NativeSolid/bin"
export PYTHONPATH=${PYTHONPATH}:${MTF_RUN}

echo ${PATH}
echo ${LD_LIBRARY_PATH}
echo ${PYTHONPATH}

