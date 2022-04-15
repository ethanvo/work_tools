#!/bin/bash

#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --sockets-per-node=2
#SBATCH --cores-per-socket=16
#SBATCH --threads-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=0-12:00:00
#SBATCH --error="%x.err"
#SBATCH --output="%x.output"
#SBATCH --account=berkelbach
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eav2136@columbia.edu

# Usage of this script:
#sbatch -J jobname job-orca-SLURM.sh , where jobname is the name of your ORCA inputfile (jobname.inp).

# Jobname below is set automatically when submitting like this: sbatch -J jobname job.sh
# Can alternatively be set manually below. job variable should be the name of the inputfile without extension (.inp)
job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})

# Activate Anaconda module
module load pyscf/gcc-openblas

# cd to scratch
cd $SLURM_SUBMIT_DIR/

# Set PYSCF max_memory
export PYSCF_MAX_MEMORY=$(bc <<< "$SLURM_MEM_PER_NODE / 1.024 / 1")

mkdir $job-tmp
# Set PySCF Temporary Directory to SLURM submit directory
export PYSCF_TMPDIR=$SLURM_SUBMIT_DIR/$job-tmp/

# Creating nodefile in scratch
echo $SLURM_NODELIST > $job.nodes

# Copy job and node info to beginning of outputfile
echo "Job execution start: $(date)" > $job.output
echo "Shared library path: $LD_LIBRARY_PATH" >>  $job.output
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  $job.output
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  $job.output
echo "Slurm node mem is: ${SLURM_MEM_PER_NODE}" >> $job.output
echo $SLURM_NODELIST >> $job.output

#Start ORCA job. ORCA is started using full pathname (necessary for parallel execution). Output file is written directly to submit directory on frontnode.
python $job.py -o $job.out
rm -rf $job-tmp/
