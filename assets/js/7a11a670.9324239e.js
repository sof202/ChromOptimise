"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[992],{8232:(e,t,i)=>{i.r(t),i.d(t,{assets:()=>l,contentTitle:()=>r,default:()=>p,frontMatter:()=>s,metadata:()=>a,toc:()=>h});var n=i(5893),o=i(1151);const s={sidebar_position:2},r="Configuration files setup",a={id:"ChromOptimise/Configuration-Files-Setup",title:"Configuration files setup",description:"You will need to create three configuration files for this pipeline to work:",source:"@site/docs/ChromOptimise/Configuration-Files-Setup.md",sourceDirName:"ChromOptimise",slug:"/ChromOptimise/Configuration-Files-Setup",permalink:"/ChromOptimise/ChromOptimise/Configuration-Files-Setup",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:2,frontMatter:{sidebar_position:2},sidebar:"documentationSidebar",previous:{title:"ChromHMM overview",permalink:"/ChromOptimise/ChromOptimise/ChromHMM-overview"},next:{title:"Pipeline explanation",permalink:"/ChromOptimise/ChromOptimise/Pipeline-Explanation"}},l={},h=[{value:"Data directory structure",id:"data-directory-structure",level:2},{value:"FilePaths.txt",id:"filepathstxt",level:2},{value:"config.R",id:"configr",level:2},{value:"LogFileManagement.sh",id:"logfilemanagementsh",level:2},{value:"ChromOptimiseConfig.txt",id:"chromoptimiseconfigtxt",level:2}];function d(e){const t={a:"a",admonition:"admonition",br:"br",code:"code",h1:"h1",h2:"h2",li:"li",p:"p",pre:"pre",ul:"ul",...(0,o.a)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(t.h1,{id:"configuration-files-setup",children:"Configuration files setup"}),"\n",(0,n.jsx)(t.p,{children:"You will need to create three configuration files for this pipeline to work:"}),"\n",(0,n.jsxs)(t.ul,{children:["\n",(0,n.jsx)(t.li,{children:"FilePaths.txt"}),"\n",(0,n.jsx)(t.li,{children:"config.R"}),"\n",(0,n.jsx)(t.li,{children:"LogFileManagement.sh."}),"\n"]}),"\n",(0,n.jsx)(t.p,{children:"These files are used by each of the R and bash scripts to aid in organisation of the scripts and avoid repetition."}),"\n",(0,n.jsx)(t.admonition,{type:"info",children:(0,n.jsx)(t.p,{children:"The scripts in this pipeline do not create the directory structure themselves. This is to avoid large files being dumped in unwanted locations. Please check your file paths are correct and then run Create_File_Structure.sh."})}),"\n",(0,n.jsx)(t.admonition,{title:"EOL errors",type:"warning",children:(0,n.jsx)(t.p,{children:"You must ensure that these files are written with EOL: \\n (LF) and not EOL: \\r\\n (CRLF)."})}),"\n",(0,n.jsxs)(t.p,{children:["After one creates each of these configuration files, place them in the 'configuration' directory. Then run the ",(0,n.jsx)(t.code,{children:"setup"})," executable from the main directory."]}),"\n",(0,n.jsx)(t.p,{children:"Note: The pipeline was completed with blueprint data in mind, if your data is already downloaded, processed, binarized etc. then the associated lines in the config files might not be required."}),"\n",(0,n.jsx)(t.h2,{id:"data-directory-structure",children:"Data directory structure"}),"\n",(0,n.jsx)(t.p,{children:"A guide for the structure of the data directory is given below (You only need to create the directories starting with an integer):"}),"\n",(0,n.jsx)(t.pre,{children:(0,n.jsx)(t.code,{className:"language-text",children:"Main_Data_Directory\n\u251c\u2500\u2500 0_Downloaded_Files\n\u251c\u2500\u2500 1_Organised_Raw_Bam_Files\n\u2502   \u251c\u2500\u2500 Epigenetic_Mark_1\n\u2502   \u251c\u2500\u2500 ...\n\u2502   \u2514\u2500\u2500 Epigenetic_Mark_n\n\u251c\u2500\u2500 2_Processed_Bam_Files\n\u2502   \u251c\u2500\u2500 Epigenetic_Mark_1\n\u2502   \u251c\u2500\u2500 ...\n\u2502   \u2514\u2500\u2500 Epigenetic_Mark_n\n\u251c\u2500\u2500 3_Subsampled_Bam_Files\n\u251c\u2500\u2500 4_Binary_Files\n\u2502   \u251c\u2500\u2500 BinSize_xxx_SampleSize_yyy\n\u2502   \u251c\u2500\u2500 ...\n\u2502   \u2514\u2500\u2500 BinSize_zzz_SampleSize_www\n\u251c\u2500\u2500 5_Model_Files\n\u251c\u2500\u2500 6_Optimum_Number_Of_States\n\u2502   \u251c\u2500\u2500 Results_From_Run_1\n\u2502   \u251c\u2500\u2500 Results_From_Run_2\n\u2502   \u251c\u2500\u2500 ...\n\u2502   \u251c\u2500\u2500 Results_From_Run_n\n\u2502   \u2514\u2500\u2500 Likelihood_Values_Of_Models\n\u251c\u2500\u2500 7_Big_Model_Files\n\u2502   \u2514\u2500\u2500 Plots\n\u2502       \u251c\u2500\u2500 Euclidean_Distance_Histrograms\n\u2502       \u2514\u2500\u2500 Transition_Maxima_Scatter_Plots\n\u2514\u2500\u2500 8_Model_Comparison_Files\n"})}),"\n",(0,n.jsx)(t.h2,{id:"filepathstxt",children:"FilePaths.txt"}),"\n",(0,n.jsx)(t.pre,{children:(0,n.jsx)(t.code,{className:"language-text",metastring:'title="FilePaths.txt"',children:'## Data directories\n\nexport MAIN_DIR="full/path/to/main/directory"\nexport DOWNLOAD_DIR="${MAIN_DIR}/path/to/downloads"\nexport RAW_DIR="${MAIN_DIR}/path/to/raw/data"\nexport PROCESSED_DIR="${MAIN_DIR}/path/to/processed/data"\nexport SUBSAMPLED_DIR="${MAIN_DIR}/path/to/subsampled/data"\nexport BINARY_DIR="${MAIN_DIR}/path/to/binary/data"\nexport MODEL_DIR="${MAIN_DIR}/path/to/chromHMM/models"\nexport OPTIMUM_STATES_DIR="${MAIN_DIR}/path/to/optimum/states/output"\nexport COMPARE_DIR="${MAIN_DIR}/path/to/comparison/files"\nexport BIG_MODELS_DIR="${MAIN_DIR}/path/to/big/models"\n\n## Script directories\n\nexport SCRIPTS_DIR="full/path/to/this/repository"\nexport RSCRIPTS_DIR="${SCRIPTS_DIR}/Rscripts"\nexport SUPPLEMENTARY_DIR="${SCRIPTS_DIR}/supplementary"\nexport LOG_DIR="${SCRIPTS_DIR}/LogFiles"\n\n## ChromHMM file locations\n\nexport CHROMHMM_MAIN_DIR="/path/to/ChromHMM/main/directory"\nexport CHROMHMM_CHROM_SIZES="${CHROMHMM_MAIN_DIR}/path/to/chromosome/sizes"\n'})}),"\n",(0,n.jsx)(t.h2,{id:"configr",children:"config.R"}),"\n",(0,n.jsxs)(t.p,{children:["To get a good value for the thresholds in the redundancy parameters section, please consult the ",(0,n.jsx)(t.a,{href:"/ChromOptimise/ChromOptimise/Supplementary-pipeline-explanation",children:"supplementary pipeline"}),"."]}),"\n",(0,n.jsx)(t.pre,{children:(0,n.jsx)(t.code,{className:"language-R",metastring:'title="config.R"',children:'## Data Directories\n\nmain_dir="path/to/main/directory"\nmodel_dir=paste0(main_dir, "path/to/model/files")\noptimum_states_dir=paste0(main_dir, "path/to/optimum/states/output")\nlikelihood_dir=paste0(optimum_states_dir, "Likelihood_Values")\ncompare_dir=paste0(main_dir, "path/to/comparison/files")\nbig_models_dir=paste0(main_dir,"path/to/big/model/files")\n\n## Plotting directories\n\ntransition_plotting_dir=paste0(big_models_dir,"/path/to/plots")\nemission_plotting_dir=paste0(big_models_dir,"/path/to/plots")\n\n## Redundancy parameters\n\nemissions_threshold=VALUE\ntransitions_threshold=VALUE\nisolation_threshold=VALUE\n\n## Number of marks used in analysis\n\nnumber_of_marks=VALUE\n'})}),"\n",(0,n.jsx)(t.h2,{id:"logfilemanagementsh",children:"LogFileManagement.sh"}),"\n",(0,n.jsxs)(t.p,{children:["This script will produce information on the time for the script to complete, reducing repetition of code\n",(0,n.jsx)(t.br,{}),"\n","(Alternatively you can use the ",(0,n.jsx)(t.code,{children:"time"})," command with the scripts, this is not available on our HPC)."]}),"\n",(0,n.jsx)(t.pre,{children:(0,n.jsx)(t.code,{className:"language-shell",metastring:'title="LogFileManagement.sh"',children:'#!/bin/bash\n## ============= ##\n##   JOB START   ##\n## ============= ##\n\necho "Job \'${SLURM_JOB_NAME}\' started at:"\ndate -u\n\nstart_time=$(date +%s)\n\nLOG_FILE_PATH="${LOG_DIR}/$SLURM_JOB_NAME/$USER"\nmkdir -p "${LOG_FILE_PATH}"\ntimestamp=$(date -u +%Y.%m.%d-%H_%M)\nexport timestamp\n\n## ============= ##\n##   FUNCTIONS   ##\n## ============= ##\n\n## ====== FUNCTION : finishing_statement() ===========================================\n## Description: Give finishing message then exit\n## Globals: \n##     start_time\n## Locals:\n##     end_time\n##     time_taken\n## Arguments:\n##     exit code\n## ===================================================================================\nfinishing_statement(){\n    echo "Job finished with exit code $1 at:"\n    date -u\n    local end_time\n    local time_taken\n    end_time=$(date +%s)\n    time_taken=$((end_time-start_time))\n    echo "Job took a total of: ${time_taken} seconds to finish."\n    exit "$1"\n}\n'})}),"\n",(0,n.jsx)(t.h2,{id:"chromoptimiseconfigtxt",children:"ChromOptimiseConfig.txt"}),"\n",(0,n.jsxs)(t.p,{children:["This is the configuration file that enables the user to run all of the files in the ",(0,n.jsx)(t.a,{href:"/ChromOptimise/ChromOptimise/Pipeline-Explanation",children:"main pipeline"})," sequentially. Options are briefly described in comments here, but for a better picture of what to put here we recommend looking at the ",(0,n.jsx)(t.a,{href:"/ChromOptimise/ChromOptimise/Pipeline-Explanation",children:"documentation"}),". If you do not plan on using certain scripts at all (which is likely the case for 0_EGADownloading.sh for example) you can just remove the options section for those selected scripts."]}),"\n",(0,n.jsx)(t.pre,{children:(0,n.jsx)(t.code,{className:"language-text",metastring:'title="ChromOptimiseConfig.txt"',children:"# Which shell script to start from (provide a number from 0 to 6)\nexport starting_script=VALUE\n\n# This is the list of marks that you intend to use in the analysis\n# Please provide this as a white space separated array\nexport LIST_OF_MARKS=(mark1 mark2 mark3 etc.)\n\n# This is a FULL file path to the FOFN which contains files you want to\n# download using the pyega3 client\nexport FILE_OF_FILE_NAMES=path/to/file\n\n# This is a threshold for the Phred score used in the processing stage\n# (which reads to discard due to low base accuracy)\nexport PRED_SCORE_THRESHOLD=VALUE\n\n# This is the sample size (as a percentage) to use in the subsampling stage\n# If your data is small in size, the recommended value is 100\nexport SAMPLE_SIZE=VALUE\n\n# This is the bin size to use during the binarization stage\n# ChromHMM recommends a default of 200\nexport BIN_SIZE=VALUE\n\n# This is the assembly that your data is alligned to.\nASSEMBLY=VALUE\n\n# This is the number of models you wish to create in the model learning stage\n# Read the documentation on 5_batch_CreateIncrementalModel.sh for help\n# here\nexport NUMBER_OF_MODELS=VALUE\n\n# This is the increment in the number of states to use between models.\n# For most cases this will likely be 1. However, if you have lots of marks\n# in your dataset a larger value might be more appropriate\nexport STATE_INCREMENT=VALUE\n"})})]})}function p(e={}){const{wrapper:t}={...(0,o.a)(),...e.components};return t?(0,n.jsx)(t,{...e,children:(0,n.jsx)(d,{...e})}):d(e)}},1151:(e,t,i)=>{i.d(t,{Z:()=>a,a:()=>r});var n=i(7294);const o={},s=n.createContext(o);function r(e){const t=n.useContext(s);return n.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function a(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(o):e.components||o:r(e.components),n.createElement(s.Provider,{value:t},e.children)}}}]);