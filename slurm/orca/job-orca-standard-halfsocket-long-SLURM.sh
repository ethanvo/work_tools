#!/bin/bash

#SBATCH --nodes=1
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=8
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --time=5-00:00:00
#SBATCH --error="%x.err"
#SBATCH --output="%x.output"
#SBATCH --account=berkelbach
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eav2136@columbia.edu

module purge
module load DefaultModules
module load orca

cd $SLURM_SUBMIT_DIR/

# Usage of this script:
#sbatch -J jobname job-orca-SLURM.sh , where jobname is the name of your ORCA inputfile (jobname.inp).

# Jobname below is set automatically when submitting like this: sbatch -J jobname job.sh
# Can alternatively be set manually below. job variable should be the name of the inputfile without extension (.inp)
job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})
# Creating nodefile in scratch

echo $SLURM_NODELIST > $SLURM_SUBMIT_DIR/$job.nodes

# Copy job and node info to beginning of outputfile

echo "Job execution start: $(date)" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Shared library path: $LD_LIBRARY_PATH" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  $SLURM_SUBMIT_DIR/$job.out
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  $SLURM_SUBMIT_DIR/$job.out
echo $SLURM_NODELIST >> $SLURM_SUBMIT_DIR/$job.out

#Start ORCA job. ORCA is started using full pathname (necessary for parallel execution). Output file is written directly to submit directory on frontnode.
$ORCA_DIR/orca $SLURM_SUBMIT_DIR/$job.inp >>  $SLURM_SUBMIT_DIR/$job.out
