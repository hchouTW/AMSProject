#!/bin/bash

# User
FL=`whoami | cut -c1`
USER=`whoami`

# VMOS
echo -e "******************************** VMOS ********************************"
voms-proxy-info --all --file ~/.ams02
echo -e "**********************************************************************"

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

proj_lst=${proj_path}/lst/${PROJFLIST}
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
if [ "${queue}" != "ams" ]; then
    echo "Queue(${queue}) is not exist."
    exit
fi

storage=${STORAGE^^}
if [ "${storage}" != "DPM" ]; then
    echo "Storage(${storage}) is not exist."
    exit
fi

storage_path=
if [ "${storage}" == "DPM" ]; then
    storage_path=${DPM_HOME}/ams02/user
    storage_check=$(dpns-ls ${storage_path} 2>&1)
    if [[ "${storage_check}" == *"No such file or directory"* ]] || [[ "${storage_check}" == *"invalid path"* ]]; then
        echo "\"${storage_path}\" is not exist."
        exit
    else
        storage_path=${storage_path}/${USER}/AMSData
        dpns-mkdir -p ${storage_path}
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


# Copy to AMSJobs
jobdir=${AMSJobs}/${PROJPATH}/${PROJVERSION}/${proj_title}
mkdir -p ${jobdir}
if [ ! -d ${jobdir} ]; then
    echo "JOBDIR(${jobdir}) is not exist."
    exit
else
    rm -rf ${jobdir}/*
fi

mkdir -p ${jobdir}/log
mkdir -p ${jobdir}/proc

cp -fa ${proj_env} ${jobdir}/env.sh
cp -fa ${proj_bin} ${jobdir}/jobexe
cp -fa ${proj_lst} ${jobdir}/flist

# libClassDef
proj_libClassDef=${proj_path}/lib/libClassDef.so
if [ -f ${proj_libClassDef} ]; then
    cp -fa ${proj_libClassDef} ${jobdir}/libClassDef.so
fi

datadir=${storage_path}/${PROJPATH}/${PROJVERSION}/${proj_title}


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
echo "#!/bin/bash
echo -e \"====  Start Run Time : \$(date)\"
echo -e \"====  Local Host : \${HOSTNAME}\"
echo -e \"====  Redhat-release  ====\"
cat /etc/redhat-release

if [ \$# -ne 3 ]; then
    echo -e \"illegal number of parameters.\"
    exit
fi

jobDir=${jobdir}
if [ ! -d \${jobDir} ]; then
    echo -e \"jobDir is not exist.\"
    exit
fi
if [ ! -d \${jobDir}/proc ]; then
    echo -e \"proc/ is not exist.\"
    exit
fi

if [ ! -d \${jobDir}/log ]; then
    echo -e \"log/ is not exist.\"
    exit
fi

source \${jobDir}/env.sh
LD_LIBRARY_PATH=\${jobDir}:\${LD_LIBRARY_PATH}

tmpData=\${jobDir}/data
mkdir -p \${tmpData}
if [ ! -d \${tmpDate} ]; then
    echo -e \"data/ is not exist.\"
    exit
fi

dataDir=${datadir}
dataDirCheck=\$(dpns-ls \${dataDir} 2>&1)
if [[ \"\${dataDirCheck}\" == *\"No such file or directory\"* ]] || [[ \"\${dataDirCheck}\" == *\"invalid path\"* ]]; then
    echo -e \"taget data/ is not exist.\"
    exit
fi


echo -e \"*************************\"
echo -e \"****  START RUNNING  ****\"
echo -e \"*************************\"

jobID=\$1
exeSatID=\$2
exeEndID=\$3

for (( exeID=\${exeSatID}; exeID<=\${exeEndID}; exeID++ ))
do
    logID=\$(printf "JOB%07i_EXE%07i" \${jobID} \${exeID})
    logLog=\${jobDir}/log/\${logID}
    logProc=\${jobDir}/proc/\${logID}
    if [ ! -f \${logProc} ]; then
        continue
    fi

    echo -e \"==== (Job) Start Time: \`date\`\\n\\n\" | tee \${logLog}
    
    ldd \${jobDir}/jobexe | tee -a \${logLog}
    
    echo -e \"\\n\\n==== (Exe) Start Time: \`date\`\" | tee -a \${logLog}
    \${jobDir}/jobexe ${event_type} \${jobDir}/flist \${exeID} ${file_per_exe} \${tmpData} 2>&1 | tee -a \${logLog}
    echo -e \"==== (Exe) End Time: \`date\`\\n\\n\" | tee -a \${logLog}
   
    FileID=\$(printf "%07i" \${exeID})
    rootFile=\`ls \${tmpData} | grep \${FileID}.root\`
    rootPath=\${tmpData}/\${rootFile}
    if [ ! -f \${rootPath} ]; then
        echo \"ROOT file is not exist.\" | tee -a \${logLog}
    else
        tagetPath=\${dataDir}/\${rootFile}
        echo -e \"==== (CopyFile) Start Time: \`date\`\" | tee -a \${logLog}
        rfcp \${rootPath} \${tagetPath}
        echo -e \"==== (CopyFile) End Time: \`date\`\\n\\n\" | tee -a \${logLog}
        tagetPathCheck=\$(dpns-ls \${tagetPath} 2>&1)
        if [[ \"\${tagetPathCheck}\" == *\"No such file or directory\"* ]] || [[ \"\${tagetPathCheck}\" == *\"invalid path\"* ]]; then
            echo "Failure." | tee -a \${logLog}
        else
            rm \${rootPath}
            rm \${logProc}
            echo "Success." | tee -a \${logLog}
        fi
    fi
    
    echo -e \"==== (Job) End Time: \`date\`\" | tee -a \${logLog}
done

echo -e \"**************************\"
echo -e \"****  FINISH RUNNING  ****\"
echo -e \"**************************\"
" > $job_script


submit_script=${jobdir}/submit.sh
echo "#!/bin/bash
echo \"**********************************************************************\"
echo \"****************************** SUBMIT ********************************\"
echo \"**********************************************************************\"

echo \"******************************** VMOS ********************************\"
voms-proxy-info --all --file ~/.ams02
echo \"**********************************************************************\"

if [ \$# -ne 1 ]; then
    echo -e \"illegal number of parameters.\"
    exit
fi

runmode=\$1
if [ \"\${runmode}\" != \"RUN\" ] && [ \"\${runmode}\" != \"RERUN\" ]; then
    echo \"\${runmode} is not correctly.\"
    exit
fi

dataDir=${datadir}
dpns-mkdir -p \${dataDir}
dataDirCheck=\$(dpns-ls \${dataDir} 2>&1)
if [[ \"\${dataDirCheck}\" == *\"No such file or directory\"* ]] || [[ \"\${dataDirCheck}\" == *\"invalid path\"* ]]; then
    echo -e \"taget data/ is not exist.\"
    exit
else
    if [ \"\${runmode}\" == \"RUN\" ]; then
        rfrm -rf \${dataDir}
        dpns-mkdir -p \${dataDir}
        dataDirCheck=\$(dpns-ls \${dataDir} 2>&1)
        if [[ \"\${dataDirCheck}\" == *\"No such file or directory\"* ]] || [[ \"\${dataDirCheck}\" == *\"invalid path\"* ]]; then
            echo -e \"taget data/ is not exist.\"
            exit
        fi
    fi
fi

jobDir=${jobdir}
if [ ! -d \${jobDir} ]; then
    echo -e \"jobDir is not exist.\"
    exit
fi

jobScript=${job_script}
if [ ! -f \${jobScript} ]; then
    echo -e \"jobScript is not exist.\"
    exit
fi

for (( jobID=0; jobID<=${totlen_job}; jobID++ ))
do
    exeSatID=\$(( ${exe_satID} + ${exe_per_job}*\${jobID} ))
    exeEndID=\$(( ${exe_satID} + ${exe_per_job}*(\${jobID}+1) - 1 ))
    if (( \${exeSatID} < ${exe_satID} )); then
        exeSatID=${exe_satID}
    fi
    if (( \${exeEndID} > ${exe_endID} )); then
        exeEndID=${exe_endID}
    fi
    
    runIt=0
    for (( exeID=\${exeSatID}; exeID<=\${exeEndID}; exeID++ ))
    do
        logID=\$(printf "JOB%07i_EXE%07i" \${jobID} \${exeID})
        if [ \${runmode} == "RUN" ]; then
            touch \${jobDir}/proc/\${logID}
            runIt=1
        else
            if [ -f \${jobDir}/proc/\${logID} ]; then
                runIt=1
            fi
        fi
    done
    
    if (( \${runIt} == 0 )); then
        continue
    fi

    jobLogID=\$(printf "JOB%07i" \${jobID})
    jobLog=\${jobDir}/log/log.\${jobLogID}

    echo \"\#!/bin/bash
sh \${jobScript} \${jobID} \${exeSatID} \${exeEndID}\" | qsub -q ${queue} -N \${jobLogID} -j oe -o \${jobLog}
done
echo \"**********************************************************************\"
echo \"************************** SUBMIT FINISH *****************************\"
echo \"**********************************************************************\"
" >> $submit_script

cfg_clear
