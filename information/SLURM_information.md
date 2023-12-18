# SLURM workload manager script elements
The scripts included in this repository were designed to be used with SLURM workload manager. As such, many elements of the scripts will not work if one does not run them through SLURM's `sbatch` command.
\
Here, we will outline the main areas of the code that depend on SLURM's `sbatch` command in order to function correctly.

## `#SBATCH` lines
At the top of every script (after the shebang) there are several comment lines that start: `#SBATCH`. These lines are read by SLURM and give meta data on how the script should be ran. If running the script normally, these lines will be interpreted as comments, resulting in loss of functionality.

## Log files
The most important aspect that is lost when running scripts outside of `sbatch` is the error logs. If one runs the scripts normally, all outputs, warnings and errors will be printed to the terminal.
\
There will also be error messages from the hard links that are created in the SET UP section of each script. Currently, the ability to create folders and useful names for these files in SLURM is not possible and limited respectively. As a result, log files are generated under temporary names and a hard link is made with a file with a more descriptive name (which allows for more organised structure).

## SLURM environment variables
Certain aspects of the scripts may lose there functionality if they depend on environment variables that are set up when running the script with SLURM's `sbatch` command. Common SLURM environment variables that were used were:
- SLURM_SUBMIT_DIR -> This is the directory that the shell was in when the script was initially ran,
- SLURM_JOB_ID -> This is the unique job identification number assigned to the job at creation,
- SLURM_JOB_NAME -> This is the name of the job as defined in the `#SBATCH` lines.

Most instances of these variables are fine as they are not required for the script to function

## Arrays
For scripts that start with 'batch', the job is designed to be ran as a SLURM array (enabling parallelisation). If one does not run these scripts as an array through `sbatch [script] --array=1-4` etc. then the script will fail. This is because the parallel processing logic that is applied in these scripts will result in division by 0 errors. Such errors occur as a blank variable is interpretted as 0 when inside of $(()).
\
The following SLURM environment variables are the reason for this behaviour:
- SLURM_ARRAY_TASK_ID
- SLURM_ARRAY_TASK_COUNT

A fair amount of effort would be required to get such scripts to work without the use of SLRUM. One would need to remove the parallelisation set up sections and ensure that all files/models are being used/manipulated correctly in subsequent sections/lines.

## Time and memory allocation
SLURM workload manager allows one to set up a peak memory allocation parameter and a maximum wall time for the submitted job. These parameters are specified in the `#SBATCH` lines at the top of the scripts. If one executes the script normally, the absence of these should not matter. One may want to implement artificial maximum wall times regardless as some programs can take a long time to finish and cause too much memory pile up (resulting in a crash).