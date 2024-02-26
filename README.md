<p align="center">
  <img src="https://github.com/sof202/ChromOptimise/assets/147140110/41a500d8-fa7d-4a3c-8270-297a16bf6b5c" />
</p>

</p>
<p align="center">
    <a href="https://img.shields.io/codefactor/grade/github/sof202/ChromOptimise" alt="CodeFactor">
        <img src="https://img.shields.io/codefactor/grade/github/sof202/ChromOptimise" /></a>
    <a href="https://img.shields.io/github/commit-activity/m/sof202/ChromOptimise" alt="Commit activity">
        <img src="https://img.shields.io/github/commit-activity/m/sof202/ChromOptimise" /></a>
    <a href="https://github.com/sof202/ChromOptimise/blob/main/LICENSE" alt="License">
        <img src="https://img.shields.io/github/license/sof202/Chromoptimise" /></a>
</p>

ChromOptimise is a pipeline that identifies the optimum number of states that should be used with [ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.)'s `LearnModel` command for a particular genomic dataset.

For more specific information, please head over to the [**wiki**](https://sof202.github.io/ChromOptimise/).

## Table of contents
- [About](#About)
- [Getting started](#getting-started)
- [Usage](#usage)
- [Software Requirements](#software-requirements)
- [Further information](#further-information)

## About
When using [ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.) to learn hidden Markov models for genomic data, it is often difficult to determine how many states to include:
- Including too many states will result in overfitting your data and introduces redundant states
- Including too few states will result in underfitting your data and thus results in lower model accuracy

This pipeline identifies the optimal number of states to use by finding a model that avoids the two above points. 

After using this pipeline, the user will have greater knowledge over their dataset in the context of ChromHMM, which will allow them to make more informed decisions as they continue to further downstream analysis.

## Getting started
1) Clone this repository
2) Ensure all [required software](#software-requirements) is installed
3) Create the configuration files using the [templates](https://sof202.github.io/ChromOptimise/ChromOptimise/Configuration-Files-Setup) provided and place them in a memorable location
4) Run the [`setup`](https://github.com/sof202/ChromOptimise/blob/main/setup) executable 
    - You may need to use `chmod +x setup` first
    - The user will be prompted for whether they want to remove lines beginning with `module` (artefact of HPC system used at UoE)
    - The user will also be prompted if they want to remove SLURM directives that are specific to the UoE HPC.

## Usage
After completing ['getting started'](#getting-started), run the master script (ChromOptimise.sh) in the command line with:
```bash
bash ChromOptimise path/to/configuration/directory
```

Alternatively, you can run each of the shell scripts in [JobSubmission](https://github.com/sof202/ChromOptimise/tree/main/JobSubmission) sequentially for each epigenetic mark. For further information and **example usage** please consult the [pipeline explanation](https://sof202.github.io/ChromOptimise/category/main-pipeline---usage-and-explanation).

Depending on your chosen dataset, you may not need to run all scripts. For example:
- If you are not downloading data from EGA, the first two scripts are not necessary 
    - Just ensure that `.bam` files are organised into directories named [[epigenetic mark name]] within the [raw files directory](https://sof202.github.io/ChromOptimise/ChromOptimise/Configuration-Files-Setup#data-directory-structure)
- If your data is already processed (quality controlled), then start from the subsampling script.
    - Again, ensure that `.bam` files are organised into directories named [[epigenetic mark name]] within the [Processed files directory](https://sof202.github.io/ChromOptimise/ChromOptimise/Configuration-Files-Setup#data-directory-structure)

There also exists supplementary scripts for further information on your chosen data set. Most importantly, thresholds used in redundancy analysis can be inferred from the results of [Redundancy_Threshold_Optimisation](https://github.com/sof202/ChromOptimise/tree/main/supplementary/Redundancy_Threshold_Optimisation). Further details for these scripts can be found in the [wiki](https://sof202.github.io/ChromOptimise/ChromOptimise/Supplementary-pipeline-explanation). 


## Software requirements
This pipeline requires a unix-flavoured OS with the following software installed:
- [Bash](https://www.gnu.org/software/bash/) (>=4.2.46(2))
- [SLURM Workload Manager](https://slurm.schedmd.com/overview.html) (>=20.02.3)
- [SAMtools](http://www.htslib.org) (>=1.9)
- [R](https://www.r-project.org) (>=3.6.0)
- [Java](https://www.java.com/en/) (>= openjdk 13.0.2)
- [ChromHMM](https://compbio.mit.edu/ChromHMM/#:~:text=ChromHMM%20is%20software%20for%20learning,and%20spatial%20patterns%20of%20marks.) (>=1.23)
- [sed](https://www.gnu.org/software/sed/) (>=4.2.2)


## Further information
This study makes use of data generated by the Blueprint Consortium. A full list of the investigators who contributed to the generation of the data is available from www.blueprint-epigenome.eu. Funding for the project was provided by the European Union's Seventh Framework Programme (FP7/2007-2013) under grant agreement no 282510 – BLUEPRINT.

For any further enquiries, please open an issue or contact Sam Fletcher:
\
s.o.fletcher@exeter.ac.uk

  
