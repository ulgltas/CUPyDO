# Set input file

export INPUT=tests/PFEM3D_Metafor/damNcomp/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/damWcomp/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/staticNcomp/input_fsi.py

# Clean output folder

rm -rf workspace
mkdir workspace

# Runs the code

export MKL_NUM_THREADS=8
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=true

python run.py -n 8 ${INPUT}