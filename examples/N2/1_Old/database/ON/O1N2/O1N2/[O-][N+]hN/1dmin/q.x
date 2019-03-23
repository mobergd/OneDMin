#!/bin/sh
#SBATCH --job-name=1dmin
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH -p bdwall
#SBATCH --time=12:00:00
#SBATCH -A CMRP
#SBATCH -o job.%j.%N.out

source ~/.bash_profile

module load molpro/2015.1_170920

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ajasper/lib/
export LD_LIBRARY_PATH
cd $SLURM_SUBMIT_DIR

mkdir -p /scratch/kmoore/1DMIN17404

cp ./ingen.pl /scratch/kmoore/1DMIN17404/.
cp -R r0 /scratch/kmoore/1DMIN17404/.

cd /scratch/kmoore/1DMIN17404

./ingen.pl 1 20
./sub-1.x
./sub-2.x

mv /scratch/kmoore/1DMIN17404 /lcrc/project/KTP/kmoore/nno/1_N2/database/ON/O1N2/O1N2/[O-][N+]hN/1dmin/.
