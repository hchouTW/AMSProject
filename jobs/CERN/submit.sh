#!/bin/bash


# Check Input Parameters
if [ $# -ne 2 ]; then
    echo "illegal number of parameters."
    exit
fi

runmode=${1^^}
if [ "${runmode}" != "RUN" ] && [ "${runmode}" != "RERUN" ]; then
    echo "${runmode} is not correctly."
    exit
fi

jobconf=$2
if [ ! -f ${jobconf} ]; then
    echo "\"${jobconf}\" is not exist."
    exit
fi

# Check Environment Variables
if [ "${AMSProj}" == "" ] || [ "${AMSCore}" == "" ] || [ "${AMSJobs}" == "" ]; then
    echo "not found (AMSProj | AMSCore | AMSJobs)"
    exit
fi

mkjob=${AMSProj}/jobs/CERN/mkjob.sh
if [ "${runmode}" == "RUN" ]; then
    sh $mkjob ${jobconf}
fi

# Ini_Parser
ini_parser=${AMSProj}/sys/shell/ini_parser.sh
if [ ! -f ${ini_parser} ]; then
    echo "\"${ini_parser}\" is not exist."
    exit
fi

source ${ini_parser}
cfg_parser $jobconf


cfgsec_PROJECT
jobdir=${AMSJobs}/${PROJPATH}/${PROJVERSION}/${PROJTITLE}
if [ ! -d ${jobdir} ]; then
    echo "JobDir(${jobdir}) is not found."
    exit
fi

submit_script=${jobdir}/submit.sh
if [ ! -f ${submit_script} ]; then
    echo "Submit Script is not found."
    exit
fi

PARAMETERS=${jobdir}/PARAMETERS
if [ ! -f ${PARAMETERS} ]; then
    echo "PARAMETERS is not found."
    exit
fi

cfgsec_QUEUE
confirm=${CONFIRM^^}
if [ "${confirm}" != "YES" ] && [ "${confirm}" != "NO" ] && [ "${confirm}" != "NONE" ]; then
    echo "Confirm(${confirm}) is not in (YES NO NONE)"
    exit
fi
echo "***************************** PARAMETERS *****************************"
cat ${PARAMETERS}
echo "******************************* RUNMODE ******************************"
echo "RUNMODE : ${runmode}" 
echo "******************************* CONFIRM ******************************"
if [ ${confirm} == "YES" ]; then
    echo "START SUBMIT JOB."
elif [ ${confirm} == "NO" ]; then
    echo "STOP SUBMIT JOB."
    exit
else
    confirm_opt=
    echo "CONFIRM (YES | NO) : "
    read confirm_opt
    case $confirm_opt in
    	"YES" | "yes" | "Y" | "y" )
    		echo "START SUBMIT JOB."
    	;;
    	"NO" | "no" | "N" | "n" )
    		echo "STOP SUBMIT JOB."
    		exit
    	;;
    	* )
    		echo "ERROR CONFIRM. EXITING ..."
    		exit
    	;;
    esac 
fi
echo "**********************************************************************"
cfg_clear

sh ${submit_script} ${runmode}
