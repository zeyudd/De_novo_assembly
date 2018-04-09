#export PATH=$PATH:/pylon5/cc3uv3p/bbrock/public/bin:/pylon5/cc3uv3p/bbrock/public/opt/bin
#export MANPATH=$MANPATH:/pylon5/cc3uv3p/bbrock/public/share/man

module use /gpfs/group/cfp102/cse597/sw/modules/
module load upc

qsub -I -A open -l nodes=1:ppn=4:scivybridge -l walltime=24:0:0 upcrun -np 3 -shared-heap=1G ./pgen test
