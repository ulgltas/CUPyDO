#!/bin/bash
# Run battery script for CUPyDO
# This is a temporary solution to remotely run the battery

# Load .profile
. ~/.profile

# Check if arguments are provided
if [ $# -ne 1 ]
  then
    echo "Usage: ./run.sh n"
    echo -e "\tn = number of threads"
    echo "ERROR: Wrong number of arguments provided!"
    exit 1
fi

# Run
echo "Running ctest on $1 threads..."
cd build
ctest -j $1

