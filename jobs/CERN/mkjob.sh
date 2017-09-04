#!/bin/bash

# User
FL=`whoami | cut -c1`
USER=`whoami`

# Check Input Parameters
if [ $# -ne 1 ]; then
    echo "illegal number of parameters."
    exit
fi
jobconf=$1
if [ ! -f ${jobconf} ]; then
    echo "\"${jobconf}\" is not exist."
    exit
fi
#runmode=${2^^}
#if [ "${runmode}" != "RUN" ] && [ "${runmode}" != "RERUN" ]; then
#    echo "${runmode} is not correctly."
#    exit
#fi


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
proj_path=${AMSCore}/${PROJPATH}/${PROJVERSION}
if [ ! -d ${proj_path} ]; then
    echo "not found \"${proj_path}\""
    exit
fi

proj_title=${PROJTITLE}
if [ "$proj_title" == "" ]; then
    proj_title=`date '+DATE%Y%m%dTIME%H%M%S'`
fi

proj_bin=${proj_path}/${PROJBIN}
if [ ! -f ${proj_bin} ]; then
    echo "not found PROJBIN(${proj_bin}\)"
    exit
fi

proj_env=${proj_path}/env/amsenv.sh
if [ ! -f ${proj_bin} ]; then
    echo "not found PROJENV(${proj_bin}\)"
    exit
fi

proj_lst=${proj_path}/job/${PROJFLIST}
if [ ! -f ${proj_lst} ]; then
    echo "not found PROJLST(${proj_lst})"
    exit
else
    lstlen=`cat ${proj_lst} | wc -l`
    if [ ${lstlen} -eq 0 ]; then
        echo "PROJLST is empty (${proj_lst})"
        exit
    fi
fi

event_type=${EVENTTYPE^^}
if [ "${event_type}" != "ISS" ] && [ "${event_type}" != "BT" ] && [ "${event_type}" != "MC" ] && [ "${event_type}" != "NONE" ]; then
    echo "TYPE(${event_type}) is not in (ISS BT MC NONE)"
    exit
fi

file_per_exe=${FILEPEREXE}
if [[ ! ${file_per_exe} =~ ^-?[0-9]+$ ]]; then
    echo "\"${file_per_exe}\" is not integer"
    exit
fi
if [ ${file_per_exe} -lt 1 ]; then
    file_per_exe=1
fi
lstlen_exe=$(( (${lstlen}/${file_per_exe}) + (${lstlen}%${file_per_exe}!=0) ))

job_region=${JOBREGION^^}
if [ "${job_region}" != "WHOLE" ] && [ "${job_region}" != "PART" ]; then
    echo "JOBREGION(${job_region}) is not in (WHOLE PART)"
    exit
fi

if [ "${job_region}" == "WHOLE" ]; then
    exe_satID=0
    exe_endID=$(( ${lstlen_exe}-1 ))
else
    exe_satID=${EXESATID}
    exe_endID=${EXEENDID}
    if [[ ! ${exe_satID} =~ ^-?[0-9]+$ ]]; then
        echo "\"${exe_satID}\" is not integer"
        exit
    fi
    if [[ ! ${exe_endID} =~ ^-?[0-9]+$ ]]; then
        echo "\"${exe_endID}\" is not integer"
        exit
    fi
    if [ ${exe_satID} -gt ${exe_endID} ]; then
        echo "satID > endID  (${exe_satID} > ${exe_endID})"
        exit
    fi
    if [ ${exe_satID} -ge ${lstlen_exe} ]; then
        echo "satID >= len  (${exe_satID} > ${lstlen_exe})"
        exit
    fi
    if [ ${exe_endID} -ge ${lstlen_exe} ]; then
        echo "endID >= len  (${exe_endID} > ${lstlen_exe})"
        exit
    fi
fi

# Input Config (QUEUE)
cfgsec_QUEUE

queue=${QUEUE}
if [ "${queue}" != "ams1nd" ]; then
    echo "Queue(${queue}) is not exist."
    exit
fi
if [ "`bqueues -u $USER | grep ${queue}`" == "" ]; then
    echo "Queue(${queue}) is not exist."
    exit
fi

storage=${STORAGE^^}
if [ "${storage}" != "EOS" ]; then
    echo "Storage(${storage}) is not exist."
    exit
fi

if [ "${storage}" == "EOS" ]; then
    storage_path=/eos/ams/user/${FL}/${USER}
    if [ ! -d ${storage_path} ]; then
        echo "\"${storage_path}\" is not exist."
        exit
    fi
fi

exe_per_job=${EXEPERJOB}
if [[ ! ${exe_per_job} =~ ^-?[0-9]+$ ]]; then
    echo "\"${exe_per_job}\" is not integer"
    exit
fi
if [ ${exe_per_job} -lt 1 ]; then
    exe_per_job=1
fi
totlen_exe=$(( ${exe_endID}-${exe_satID}+1 ))
totlen_job=$(( ${totlen_exe}/${exe_per_job} + (${totlen_exe}%${exe_per_job}!=0) ))

confirm=${CONFIRM^^}
if [ "${confirm}" != "YES" ] && [ "${confirm}" != "NO" ] && [ "${confirm}" != "NONE" ]; then
    echo "Confirm(${confirm}) is not in (YES NO NONE)"
    exit
fi

echo "**********************************************************************
***************************** PARAMETERS *****************************
**********************************************************************
[PROJECT]
PROJPATH        ${PROJPATH}
PROJVERSION     ${PROJVERSION}
PROJTITLE       ${proj_title}
PROJBIN         ${PROJBIN}
PROJFLST        ${PROJFLIST}
EVENTTYPE       ${event_type}
FILEPEREXE      ${file_per_exe}
JOBREGION       ${job_region}

[QUEUE]
QUEUE           ${queue}
STORAGE         ${storage}
EXEPERJOB       ${exe_per_job}

[JOBS]
EXESATID        ${exe_satID}
EXEENDID        ${exe_endID}
TOTALEXE        ${totlen_exe}
TOTALJOB        ${totlen_job}
**********************************************************************"
echo "******************************** CONFIRM *******************************"
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


# Copy to AMSJobs
jobdir=${AMSJobs}/${PROJPATH}/${PROJVERSION}/${proj_title}
mkdir -p ${jobdir}
if [ ! -d ${jobdir} ]; then
    echo "JOBDIR(${jobdir}) is not exist."
    exit
fi

mkdir -p ${jobdir}/log
mkdir -p ${jobdir}/proc

cp -fa ${proj_env} ${jobdir}/env.sh
cp -fa ${proj_bin} ${jobdir}/jobexe
cp -fa ${proj_lst} ${jobdir}/flist

echo "[PROJECT]
PROJPATH        ${PROJPATH}
PROJVERSION     ${PROJVERSION}
PROJTITLE       ${proj_title}
PROJBIN         ${PROJBIN}
PROJFLST        ${PROJFLIST}
EVENTTYPE       ${event_type}
FILEPEREXE      ${file_per_exe}
JOBREGION       ${job_region}
EXESATID        ${exe_satID}
EXEENDID        ${exe_endID}

[QUEUE]
QUEUE           ${queue}
STORAGE         ${storage}
EXEPERJOB       ${exe_per_job}

[JOBS]
TOTALEXE        ${totlen_exe}
TOTALJOB        ${totlen_job}" > ${jobdir}/PARAMETERS

# Job Script
job_script=${jobdir}/job.sh
echo "#!bin/bash
#shopt -s -o nounset
EOScomd='/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
echo -e \"====  Start Run Time : '$(date)'  ====\"
echo -e \"====  Local Host : "'$HOSTNAME'"  ====\"
echo -e \"====  Redhat-release  ====\"
cat /etc/redhat-release

echo -e \"*************************\"
echo -e \"****  START RUNNING  ****\"
echo -e \"*************************\"

JobID=$1
ExeSatID=$2
ExeEndID=$3

for (( exeID=${ExeSatID}; exeID<=${ExeEndID}; exeID++ ))
do

    # check JOB{ID}EXE{ID} exist


    echo -e \"****  RUNNING  ****\" | tee ${logID_File}
    ./jobexe #{EventMode} ${Stream} ${exeID} #{FilePerExe} ${TmpVersionStreamSrgDir} 2>&1 | tee -a ${logID_File}

    rootFile=`ls ${TmpVersionStreamSrgDir} | grep ${fileName}`
    cp ${TmpVersionStreamSrgDir}/$rootFile ${VersionStreamSrgDir}/$rootFile

	/bin/rm ${TmpVersionStreamSrgDir}/$rootFile
	echo -e \"****  FINISH RUNNING  ****\" | tee -a ${logID_File}
done

echo -e \"**************************\"
echo -e \"****  FINISH RUNNING  ****\"
echo -e \"**************************\"
" > $job_script


#submit_script=${jobdir}/submit.sh
#echo "#!bin/bash
#runmode=$1
#
##if run; then touch proc/JOB{ID}EXE{ID}
##cp script exe flist to host
##cp proc/JOB{ID}EXE{ID} to host (if rurun; then cp exist JOB{ID}EXE{ID} to host)
#if rerun; check files.....
#
#
#" >> $submit_script





#TOTAL_RUN_EXE      #{TotalExe}
#EXE_START_ID       #{ExeStartID}
#EXE_END_ID         #{ExeEndID}
#
#TOTAL_RUN_JOB      #{TotalJob}
#JOB_START_ID       #{JobStartID}
#JOB_END_ID         #{JobEndID}
