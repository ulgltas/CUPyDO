# Set input file

# export INPUT=tests/PFEM3D_Metafor/birdAxisym/input_fsi.py
export INPUT=tests/PFEM3D_Metafor/staticNcomp/input_fsi.py
export INPUT=tests/PFEM3D_Metafor/staticNcompConsist/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/damNcomp/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/damWcomp/input_fsi.py

# export INPUT=tests/PFEM_Metafor/beam/beam_fsi.py
# export INPUT=tests/PFEM_Metafor/bird/impactaxi_fsi.py
# export INPUT=tests/PFEM_Metafor/bird/impact_fsi.py
# export INPUT=tests/PFEM_Metafor/static/cylinder_fsi.py
# export INPUT=tests/PFEM_Metafor/bird/lsdyna_fsi.py
# export INPUT=tests/PFEM_Metafor/wcolumn/elasticgate/rho1100_IQNILS6_fsi.py
# export INPUT=tests/PFEM_Metafor/wcolumn/elasticgate/rho1100_IQNILS_fsi.py
# export INPUT=tests/PFEM_Metafor/wcolumn/elasticgate/rho1100_fsi.py
# export INPUT=tests/PFEM_Metafor/wcolumn/obstacle/matching_fsi.py
# export INPUT=tests/PFEM_Metafor/wcolumn/obstacle/nonmatching_fsi.py

# Clean output folder

rm -rf workspace
mkdir workspace

# Runs the code

export MKL_NUM_THREADS=8
export OMP_NUM_THREADS=8
python3 run.py -k 8 ${INPUT}
