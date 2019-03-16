#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition BioCompute
#SBATCH --nodes=1
#SBATCH --ntasks=1  # used for MPI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=12  # cores per task
#SBATCH --mem-per-cpu=8G  # memory per core (default is 1GB/core)
#SBATCH --time 0-02:00  # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=kinglab  # investors will replace this with their account name
# once we get you added this ^ line should be '--account=kinglab'
## labels and outputs
#SBATCH --job-name=LearnAbundance_eb
#SBATCH --output=test_%A_%a.out  # Standard output
#SBATCH -e error_%A_%a.err # Standard error

## notifications
#SBATCH --mail-user=pawtk4@mail.missouri.edu  # email address for notifications
#SBATCH --mail-type=END,FAIL  # which type of notifications to send
#-------------------------------------------------------------------------------
module load stringtie/stringtie-1.3.3b 
 
 
echo "### Starting at: $(date) ###"
 
# load modules then display what we have
module load  

 
# Science goes here:

COMMANDA=`head -n ${SLURM_ARRAY_TASK_ID} ../rna-seq_data_bam/S08_Abundances_LearnMemRNA_eb.txt | tail -n 1` 
$COMMANDA
 
echo "### Ending at: $(date) ###"

