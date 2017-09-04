#!/bin/bash

# User
FL=`whoami | cut -c1`
USER=`whoami`

# Check Input Parameters
if [ $# -lt 1 ]; then
    echo "illegal number of parameters."
    exit
fi
jobconf=$1
if [ ! -f ${jobconf} ]; then
    echo "\"${jobconf}\" is not exist."
    exit
fi
runmode=${2^^}
if [ "${runmode}" != "RUN" ] && [ "${runmode}" != "RERUN" ]; then
    echo "${runmode} is not correctly."
    exit
fi

# Check Environment Variables
if [ "${AMSProj}" == "" ] || [ "${AMSCore}" == "" ] || [ "${AMSJobs}" == "" ]; then
    echo "not found (AMSProj | AMSCore | AMSJobs)"
    exit
fi


# Ini_Parser
ini_parser=${AMSProj}/sys/shell/ini_parser.sh
if [ ! -f ${ini_parser} ]; then
    echo "\"${ini_parser}\" is not exist."
    exit
fi

source ${ini_parser}
cfg_parser $jobconf

# Input Config (PROJECT)
cfgsec_PROJECT
submit_script=${AMSCore}/${PROJPATH}/${PROJVERSION}/submit.sh

sh ${submit_script} ${runmode}




confirm=${CONFIRM^^}
if [ "${confirm}" != "YES" ] && [ "${confirm}" != "NO" ] && [ "${confirm}" != "NONE" ]; then
    echo "Confirm(${confirm}) is not in (YES NO NONE)"
    exit
fi
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
