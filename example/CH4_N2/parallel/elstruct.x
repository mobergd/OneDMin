module load molpro/2015.1_170920
molpro -n 4 --nouse-logfile --no-xml-output -L /soft/molpro/2015.1_170920/bebop/molprop_2015_1_linux_x86_64_i8/lib/ -d /scratch/$USER -o qc.out -s qc.in
