#!/bin/bash

# User
FL=`whoami | cut -c1`
USER=`whoami`

# Check Environment
if [ $# -ne 2 ]; then
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

if [ "${AMSProj}" == "" ] || [ "${AMSCore}" == "" ] || [ "${AMSJobs}" == "" ]; then
    echo "not found (AMSProj | AMSCore | AMSJobs)"
    exit
fi

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

exe_per_job=${EXEPERJOB}
if [[ ! ${} =~ ^-?[0-9]+$ ]]; then
    echo "\"${exe_per_job}\" is not integer"
    exit
fi
if [ ${exe_per_job} -lt 1 ]; then
    exe_per_job=1
fi
totlen_exe=$(( ${exe_endID}-${exe_satID}+1 ))
totlen_job=$(( ${totlen_exe}/${exe_per_job} + (${totlen_exe}%${exe_per_job}!=0) ))

storage=${STORAGE^^}
if [ "${queue}" != "EOS" ]; then
    echo "Storage(${storage}) is not exist."
    exit
fi

if [ "${queue}" == "EOS" ]; then
    storage_path=/eos/ams/user/${FL}/${USER}
    if [ ! -d ${storage_path} ]; then
        echo "\"${storage_path}\" is not exist."
        exit
    fi
fi

confirm=${CONFIRM^^}
if [ "${confirm}" != "YES" ] || [ "${confirm}" != "NO" ] || [ "${confirm}" != "NONE" ]; then
    echo "Confirm(${confirm}) is in (YES NO NONE)"
    exit
fi

echo "========== PARAMETERS ==========
[PROJECT]
PATH        ${PROJPATH}/${PROJVERSION}
TITLE       ${proj_title}
BIN         ${proj_bin}
FLST        ${proj_lst}
EVENTTYPE   ${event_type}
REGION      ${job_region}
EXESATID    ${exe_satID}
EXEENDID    ${exe_endID}

[QUEUE]
QUEUE       ${queue}
EXEPERJOB   ${exe_per_job}
STORAGE     ${storage}
======================================="
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
echo "================================================="



jobdir=${AMSJobs}/${PROJPATH}/${PROJVERSION}/${proj_title}
mkdir -p ${jobdir}
if [ ! -d ${jobdir} ]; then
    echo "JOBDIR(${jobdir}) is not exist."
    exit
fi

mkdir -p ${jobdir}/log
mkdir -p ${jobdir}/proc

cp -fa ${proj_env} ${jobdir}/env.sh
cp -fa ${proj_path}/${proj_bin} ${jobdir}/jobexe
cp -fa ${proj_path}/${proj_lst} ${jobdir}/flist

echo "PARAMETERS" > ${jobdir}/PARAMETERS

job_script=${jobdir}/job.sh
echo "
" > $job_script




#ParamFileName = "PARAMETERS"
#ParamFilePath = "#{VersionStreamLogDir}/#{ParamFileName}"
#ParamFile = File.open("#{ParamFilePath}", "w")
#ParamFile << """
#QUEUE_MODE         #{QueueMode}
#STORAGE_MODE       #{StorageMode}
#PROJECT_MODE       #{ProjectMode}
#"""
#if (ProjectMode.eql? "ANALYSIS"); then
#ParamFile << """
#SUBPROJECT_MODE    #{SubProjectMode}
#"""
#end
#ParamFile << """
#RUN_MODE           #{RunMode}
#EVENT_MODE         #{EventMode}
#
#VERSION            #{Version}
#VERSION_TITLE      #{VersionTitle}
#VERSION_STREAM     #{VersionStream}
#JOB_EXE            #{JobExe}
#STREAM             #{Stream}
#
#FILE_PER_EXE       #{FilePerExe}
#EXE_PER_JOB        #{ExePerJob}
#
#TOTAL_RUN_EXE      #{TotalExe}
#EXE_START_ID       #{ExeStartID}
#EXE_END_ID         #{ExeEndID}
#
#TOTAL_RUN_JOB      #{TotalJob}
#JOB_START_ID       #{JobStartID}
#JOB_END_ID         #{JobEndID}
#
#"""
#ParamFile.close
