# Set input file

# export INPUT=tests/PFEM3D_Metafor/birdAxisym/input_fsi.py
export INPUT=tests/PFEM3D_Metafor/clampedBeam/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/clampedBeamConsist/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/damBreak/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/damBreakWeakly/input_fsi.py
# export INPUT=tests/PFEM3D_Metafor/pureConduction/input_fsi.py

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
