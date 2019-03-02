#!/bin/bash

cwd=`pwd`
scd="/scratch/$USER/1DMIN$$"

if [ "$nam" == "" ]; then
        nam='job'
fi

cat << EOF > q.x
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

LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/home/ajasper/lib/
export LD_LIBRARY_PATH
cd \$SLURM_SUBMIT_DIR

mkdir -p $scd

cp ./ingen.pl $scd/.
cp -R r0 $scd/.

cd $scd

./ingen.pl 1 20
./sub-1.x
./sub-2.x

mv $scd $cwd/.
EOF

chmod u+x q.x

sbatch q.x
