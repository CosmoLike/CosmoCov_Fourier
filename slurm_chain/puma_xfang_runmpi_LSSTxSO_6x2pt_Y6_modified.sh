#!/bin/bash
#SBATCH --job-name=lsstso6_6x2
#SBATCH --account=timeifler
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --nodes=12
#SBATCH --ntasks=1128
#SBATCH --ntasks-per-node=94
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=25:00:00

#SBATCH --output=%A.out
#SBATCH --error=%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xfang@email.arizona.edu

cd $SLURM_SUBMIT_DIR

module load gsl/2.6
module load python/3.6/3.6.5
module load openmpi3/3.1.4

### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
export LD_LIBRARY_PATH="/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib:$LD_LIBRARY_PATH"
date
/usr/bin/time mpirun -x LD_LIBRARY_PATH --mca pml ob1 --mca btl ^openib -n 1128 python3 run_sampler/xfang_runLSSTxSO_Y6_nx2pt_modified.py 2 1128
date