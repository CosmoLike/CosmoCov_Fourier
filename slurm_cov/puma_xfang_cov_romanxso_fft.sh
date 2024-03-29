#!/bin/bash
#SBATCH --job-name=RxSO_cv
#SBATCH --account=timeifler
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --nodes=1
#SBATCH --ntasks=15
#SBATCH --array=0-658
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=15:00:00
#SBATCH --output=%A.out
#SBATCH --error=%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=xfang@email.arizona.edu

cd $SLURM_SUBMIT_DIR
module load gsl/2.6
module load python/3.6/3.6.5
N=9870
for ((i = 0; i < SLURM_NTASKS; i++)); do
	i_cov=$((SLURM_ARRAY_TASK_ID*SLURM_NTASKS+i+1));
	if [[ ${i_cov} -gt N ]]; then break; fi;
	srun --ntasks 1 --exclusive -c 1 ./compute_covariances_fourier_10x2pt_fft ${i_cov} 3 >&/home/u1/xfang/output/jb_covroman_${i_cov}.log
done;