eval "$(conda shell.bash hook)"
conda activate
python equi_ellip.py
./clean.sh
./ibmc
