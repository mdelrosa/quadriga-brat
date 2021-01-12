#!/bin/bash -l

#SBATCH --mail-type=ALL
#SBATCH --mail-user=mdelrosa@ucdavis.edu
#SBATCH -p med
#SBATCH --time=4-0 # 4 days-0 hours
#SBATCH -c 1  # 1 cpu cores per node
#SBATCH --mem=16000 # use max node
# set a=batch_glob, then submit with the following line to get properly named jobs:
# sbatch --export=a=$a,gen_dir=$b --job-name=outslow_$a data_gen_outdoor_slow.sh

module load matlab/2019a
matlab -r "addpath('${quadriga_src}'); cd ${gen_dir}; batch_num=${batch_num}; C01_3GPP_36873_3D_UMa_NLOS();"