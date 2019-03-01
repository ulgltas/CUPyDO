#!/bin/bash
# Run script for CUPyDO
# External solver dir should place next to CUPyDO dir

# Load .profile
. ~/.profile

# Set paths
CUP_DIR="$(dirname "$0")"
TOP_DIR=$(cd "${CUP_DIR}/.."; pwd)
MTF_RUN="${TOP_DIR}/Metafor"
# LIB
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MTF_RUN}
# PYTHON
export PYTHONPATH=${PYTHONPATH}:${CUP_DIR}
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/SU2/bin"
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/waves"
export PYTHONPATH="${PYTHONPATH}:${TOP_DIR}/NativeSolid/bin"
export PYTHONPATH=${PYTHONPATH}:${MTF_RUN}

# Print paths
echo PATH=${PATH}
echo LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
echo PYTONPATH=${PYTHONPATH}
echo ""

# Check if arguments are provided
if [ $# -ne 2 ]
  then
    echo "Usage: ./run.sh f n"
    echo -e "\tf = path/to/file (run single test) || f = all (run ctest)"
    echo -e "\tn = number of threads"
    echo "ERROR: Wrong number of arguments provided!"
    exit 1
fi

# Run
if [ "$1" == "all" ]
  then
    echo "Running ctest on $2 threads..."
    cd build
    ctest
else
    if [[ -f $1 ]]
      then
        echo "Running file $1 on $2 threads..."
        python $1 -n $2
    else
        echo "$1 is not a file"
        exit 1
    fi
fi

