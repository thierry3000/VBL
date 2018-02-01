#!/bin/bash

#SBATCH --job-name=testSLURM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=2GB
#SBATCH --nodelist=leak61
echo #SBATCH --exclude=leak[57-64]

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

#cd /localdisk/thierry/buildopt_vbl/tools/tests/slurm
srun /localdisk/thierry/buildopt_vbl/tools/tests/slurm/test_omp_slurm

#end of run_script.sh
