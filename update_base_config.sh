#!/bin/bash

PAM_BASE_CONFIG_URL='https://raw.githubusercontent.com/sanger-pathogens/nextflow-commons/master/configs/nextflow.config'
PAM_BASE_CONFIG_OUTFILE='pam_generic_base.config'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PAM_BASE_CONFIG_OUTFILE_PATH=${SCRIPT_DIR}/configs/${PAM_BASE_CONFIG_OUTFILE}

wget -O ${PAM_BASE_CONFIG_OUTFILE_PATH} ${PAM_BASE_CONFIG_URL}

sed -i 's/standard {/sanger_lsf {/g' $PAM_BASE_CONFIG_OUTFILE_PATH
