# Analyze
module load 2025
module load GROMACS/2026.3-foss-2025a-CUDA-12.8.0
sh analysis/prepare.py
python analysis/rmsd.py
python analysis/rmsf.py
python analysis/rg.py
python analysis/pca.py