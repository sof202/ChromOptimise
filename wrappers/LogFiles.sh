#!/bin/bash

## ===================== ##
##   LOG FILE MOVEMENT   ##
## ===================== ##

LOG_FILE_PATH="${LOG_DIR}/${USER}/${SLURM_JOB_NAME}"
mkdir -p "${LOG_FILE_PATH}"
timestamp=$(date +%d-%h~%H-%M)

if [[ -n ${SLURM_ARRAY_TASK_ID} ]]; then
    mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
    "${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.log"
    mv "${SLURM_SUBMIT_DIR}/temp${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
    "${LOG_FILE_PATH}/${SLURM_ARRAY_JOB_ID}~${SLURM_ARRAY_TASK_ID}~$timestamp.err"
else
    mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.log" \
    "${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.log"
    mv "${SLURM_SUBMIT_DIR}/temp${SLURM_JOB_ID}.err" \
    "${LOG_FILE_PATH}/${SLURM_JOB_ID}~$timestamp.err"
fi

## ================= ##
##   START MESSAGE   ##
## ================= ##

cat << START_MESSAGE
Job ${SLURM_JOB_NAME} started at:
$(date -u)

File paths used:
$(cat "${configuration_directory}/FilePaths.txt")

R config used:
$(cat "${configuration_directory}/config.R")

START_MESSAGE

start_time=$(date +%s)
export start_time

## =============== ##
##   END MESSAGE   ##
## =============== ##

## ====== FUNCTION : finishing_statement() ===========================================
## Description: Give finishing message then exit
## Globals: 
##     start_time
## Locals:
##     end_time
##     time_taken
## Arguments:
##     exit code
## ===================================================================================
finishing_statement(){
    local end_time
    local time_taken
    end_time=$(date +%s)
    time_taken=$((end_time-start_time))
cat << END_MESSAGE
Job finished with exit code $1 at:
date -u

Job took a total of: ${time_taken} seconds to finish.
END_MESSAGE
    exit "$1"
}
