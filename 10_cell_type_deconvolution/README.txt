Made based on guide from https://github.com/smangul1/rop/wiki/ROP-analysis-via-qsub, but with substantial changes to the commands to accomondate project structure
Files/data related to the running of the pipeline and not actual results are stored in the scripts folder rather than the data folder.

############################
# STEP 1: RUN RUP PIPELINE #
############################

To run the pipeline, the following modules are needed (These are loaded with each script/job):
anaconda3/2023.09-0
samtools/1.8
ncbi-blast/2.15.0+
bwa/0.7.17

The files for running are made like so:

### BATCH 1 ###
# Make file containing all file names needed for analyses
ls ../../data/08_unmapped_read_mapping_to_viral_genomes/batch1_pass2_unmapped/*Unmapped.out.mate* | awk -F "batch1_pass2_unmapped/" '{print $2}'>sample.txt

# Make folders for run files and results
mkdir batch1_run_files
mkdir ../../data/08_unmapped_read_mapping_to_viral_genomes/results_rop_batch1

# Make run files for all fastqs
while read line; do echo -e "module load tools\nmodule load anaconda3/2023.09-0\nmodule load samtools/1.8\nmodule load perl\nmodule load ncbi-blast/2.15.0+\nmodule load bwa/0.7.17\ncd /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/scripts/08_unmapped_read_mapping_to_viral_genomes/batch1_run_files \nbash /home/projects/ssi_gen1/data/pregnancy/unmapped/rop_latest/rop/rop_frehed_v2.sh -s all -f -i /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/data/08_unmapped_read_mapping_to_viral_genomes/batch1_pass2_unmapped/${line} /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/data/08_unmapped_read_mapping_to_viral_genomes/results_rop_batch1/${line}" > batch1_run_files/run_${line}.sh; done<samples_batch1.txt

# Make file with all the qsub commands
ls batch1_run_files/run*sh | awk '{run_file=$1; sub(/^.*\/run_/, "", $1); sub(/\.sh$/, "", $1); print "qsub -W group_list=ssi_gen1 -A ssi_gen1 -N ROP_"$1" -l nodes=1:thinnode:ppn=10,mem=16GB,walltime=00:30:00 "run_file}' > qsub_commands_all_batch1.sh

# All jobs are then submitted to the queue like so:
nohup /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/scripts/08_unmapped_read_mapping_to_viral_genomes/qsub_commands_all_batch1.sh

NOTE: If you want to queue a smaller number of files at a time, you can split the qsub_commands_all_batch1.sh into smaller segments (here shown 250) like so:
split -l 250 qsub_commands_all_batch1.sh qsub_commands_all_batch1_ 


### BATCH 2 ###
# Make file containing all file names needed for analyses
ls ../../data/08_unmapped_read_mapping_to_viral_genomes/batch2_pass2_unmapped/*Unmapped.out.mate* | awk -F "batch2_pass2_unmapped/" '{print $2}'>samples_batch2.txt

# Make folders for run files and results
mkdir batch2_run_files
mkdir ../../data/08_unmapped_read_mapping_to_viral_genomes/results_rop_batch2

# Make run files for all fastqs
while read line; do echo -e "module load tools\nmodule load anaconda3/2023.09-0\nmodule load samtools/1.8\nmodule load perl\nmodule load ncbi-blast/2.15.0+\nmodule load bwa/0.7.17\ncd /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/scripts/08_unmapped_read_mapping_to_viral_genomes/batch2_run_files \nbash /home/projects/ssi_gen1/data/pregnancy/unmapped/rop_latest/rop/rop_frehed_v2.sh -s all -f -i /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/data/08_unmapped_read_mapping_to_viral_genomes/batch2_pass2_unmapped/${line} /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/data/08_unmapped_read_mapping_to_viral_genomes/results_rop_batch2/${line}" > batch2_run_files/run_${line}.sh; done<samples_batch2.txt

# Make file with all the qsub commands
ls batch2_run_files/run*sh | awk '{run_file=$1; sub(/^.*\/run_/, "", $1); sub(/\.sh$/, "", $1); print "qsub -W group_list=ssi_gen1 -A ssi_gen1 -N ROP_"$1" -l nodes=1:thinnode:ppn=10,mem=16GB,walltime=01:00:00 "run_file}' > qsub_commands_all_batch2.sh

# All jobs are then submitted to the queue like so:
nohup /home/projects/ssi_gen1/data/pregnancy/finalProjectAnalyses/scripts/08_unmapped_read_mapping_to_viral_genomes/qsub_commands_all_batch2.sh


################################
# STEP 2: EXTRACT VIRAL COUNTS #
################################
For each batch, the script to extract a viral counts file is made:
08b_viralCounts_batch*.sh

The rest can then be run according to script naming.


