#!/bin/sh
######## Slurm header #########
#SBATCH --job-name=plots_MFCblockInfo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --output="%x.o%"
#SBATCH --mail-type=BEGIN,END,FAIL

##### print information #####
echo "Working Directory: 	$PWD"
echo "Running on host		`hostname`"
echo "Job id:				$SLURM_JOB_ID"
echo "Job name:				$SLURM_JOB_NAME"
echo "Number of cores allocated to job: $SLURM_NNODES"
echo "Number of nodes allocated to job: $SLURM_NTASKS"

# load modules
module load compiler/gnu/10.2
module load math/R_gnu/3.6.3-gnu-10.2

# Run R script
R CMD BATCH --no-save --no-restore --slave plots_simulated.R
