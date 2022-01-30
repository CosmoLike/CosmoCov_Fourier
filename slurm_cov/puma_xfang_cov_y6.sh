#!/bin/bash
#SBATCH --job-name=LSSTxSO_cov
#SBATCH --account=timeifler
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --nodes=11
#SBATCH --ntasks=1000
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=10:00:00
#SBATCH --output=%A.out
#SBATCH --error=%A.err
cd $SLURM_SUBMIT_DIR
module load gsl/2.6
module load python/3.6/3.6.5
N=11781
for i in {0..999}; do
	i_cov=$((SLURM_ARRAY_TASK_ID*1000+i+1));
	if [[ ${i_cov} -gt N ]]; then break; fi;
	srun --nodes 1 --ntasks 1 ./compute_covariances_fourier ${i_cov} 1 >&/home/u1/xfang/output/job_output_${i_cov}.log
done;