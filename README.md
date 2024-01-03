# Optimum States Analysis

## Finds the optimum number of states that should be used with ChromHMM's `LearnModel` with a particular data set

### About
This repository contains bash scripts and R scripts that aim to identify an 'intelligent' number of states to include in the hidden Markov models produced by ChromHMM's `LearnModel` command (as opposed to guessing this value).
\
When producing models with ChromHMM, it is often difficult to detemine how many states to include:
- Including too many states will result in overfitting your data and introduces redundant states
- Including too few states will result in underfitting your data and thus loses model accuracy

This pipeline produces many different models using ChromHMM and aims to determine the optimal number of states by finding a model that avoids the two above points. 
\
\
The optimum number of states to use will likely change depending on different data sets and user inputs. Examples of factors that can change the optimum number of states for the inputed data set are:
- The Phred score threshold 
    - Which reads are rejected during alignment process due to low sequencing accuracy
- The sample size
    - The pipeline includes a subsampling stage that can help in identifying the completeness of your data
- The bin size used in the binarization process (binarization is a part of the ChromHMM pipeline)
    - Smaller bin sizes will be more suceptible to random background noise in the data
    - Larger bin sizes lose data precision
- The number and types of epigenetic marks included in the data set
- The random noise apparent in the raw sequencing data
- The presence/absence of a control file

After using this pipeline, the user will have greater knowledge over their dataset in the context of ChromHMM, which will allow them to make more informed decisions as they continue to further downstream analysis.
\
\
Note: This pipeline was designed with a very large dataset in mind. The dataset being blueprint from EGA (which includes a large amount of ChIP-Seq data from Mature Neutrophils in blood samples). Keep this in mind when experimenting with subsampling.
\
\
Note 2: For a quick explanation on how a script works, run the script with the `-h` or `--help` flags.
### Required software
This pipeline requires the following software:
- Bash (>=4.2.46(2))
- SLURM Workload Manager (>=20.02.3)
- SAMtools (>=1.9)
- R (>=4.3.1)
- Java openjdk (>=13.0.2)
- ChromHMM (>=1.23)

### Further information
Further information on the following topics can be found in the [wiki](https://github.com/sof202/ChromHMM_Optimum_States/wiki)
- Pipeline_Explanation.md -> Script / pipeline specific information
- Optimal_States_Analysis_Pipeline.pdf -> Pipeline process in schematic representation
- Processing_Times.md -> Estimated processing times for scripts subject to CPU performance 
- Memory_Profiling.md -> Estimated peak heap memory consumption for scripts
- SLURM_Information.md -> Information on the SLURM workload manager used for job submission
- Config_Setup.md -> Provides templates for the config files sourced in the scripts.

For any further enquiries, please contact:
\
Sam Fletcher: s.o.fletcher@exeter.ac.uk

  