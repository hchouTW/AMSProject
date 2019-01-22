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

# +JobFlavour = "tomorrow"
# espresso     20m
# microcentury 1h
# longlunch    2h
# workday      8h
# tomorrow     1d
# testmatch    3d
# nextweek     7d
queue=${QUEUE}
if [ "${queue}" != "espresso" ] && [ "${queue}" != "microcentury" ] && [ "${queue}" != "longlunch" ] && [ "${queue}" != "workday" ] && [ "${queue}" != "tomorrow" ] && [ "${queue}" != "testmatch" ] && [ "${queue}" != "nextweek" ]; then
    echo "Queue(${queue}) is not exist."
    exit
fi


storage=${STORAGE^^}
if [ "${storage}" != "EOS" ]; then
    echo "Storage(${storage}) is not exist."
    exit
fi

storage_path=
if [ "${storage}" == "EOS" ]; then
    storage_path=/eos/ams/user/${FL}/${USER}
    if [ ! -d ${storage_path} ]; then
        echo "\"${storage_path}\" is not exist."
        exit
    else
        storage_path=${storage_path}/AMSData
        mkdir -p ${storage_path}
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

mkdir -p ${jobdir}/out
mkdir -p ${jobdir}/err
mkdir -p ${jobdir}/log
if [ ! -f ${jobdir}/proc ]; then
    touch ${jobdir}/proc
fi

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
EOScomd='/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select'
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

procFile=\${jobDir}/proc
if [ ! -f \${procFile} ]; then
    echo -e \"proc is not exist.\"
    exit
fi

if [ ! -d \${jobDir}/out ]; then
    echo -e \"out/ is not exist.\"
    exit
fi

if [ ! -d \${jobDir}/err ]; then
    echo -e \"err/ is not exist.\"
    exit
fi

if [ ! -d \${jobDir}/log ]; then
    echo -e \"log/ is not exist.\"
    exit
fi

cp -fa \${jobDir}/env.sh \${PWD}/env.sh
cp -fa \${jobDir}/jobexe \${PWD}/jobexe
cp -fa \${jobDir}/flist  \${PWD}/flist

# libClassDef
libClassDef=\${jobDir}/libClassDef.so
if [ -f \${libClassDef} ]; then
    cp -fa \${libClassDef} \${PWD}/libClassDef.so
fi

source \${PWD}/env.sh
LD_LIBRARY_PATH=\${PWD}:\${LD_LIBRARY_PATH}

source /eos/ams/user/h/hchou/ExternalLibs/BIN/CERN_TRACKSys.sh

tmpData=\${PWD}/data
mkdir -p \${tmpData}
if [ ! -d \${tmpDate} ]; then
    echo -e \"data/ is not exist.\"
    exit
fi

dataDir=${datadir}
if [ ! -d \${dateDir} ]; then
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
    logID=\$(printf \"JOB%07i_EXE%07i\" \${jobID} \${exeID})
    if [ \"\`grep -n \${logID} \${procFile}\`\" != \"\" ]; then 
        continue
    fi

    echo -e \"==== (Job) Start Time: \`date\`\\n\\n\"
    
    echo \"PWD:  \${PWD}\"
    ls -alh \${PWD}
    echo -e \"\\n\\n\"
    ldd \${PWD}/jobexe
    
    echo -e \"\\n\\n==== (Exe) Start Time: \`date\`\"
    ./jobexe ${event_type} flist \${exeID} ${file_per_exe} \${tmpData}
    echo -e \"==== (Exe) End Time: \`date\`\\n\\n\"
    
    FileID=\$(printf \"%i\" \${exeID})
    rootFile=\`ls \${tmpData} | grep \${FileID}.root\`
    rootPath=\${tmpData}/\${rootFile}
    
    echo -e \"\\n\\n\"
    ls -alh \${tmpData}
    echo -e \"\\n\"
    ls -alh \${rootPath}
    echo -e \"\\n\\n\"
    
    if [ ! -f \${rootPath} ]; then
        echo \"ROOT file is not exist.\"
    else
        tagetPath=\${dataDir}/\${rootFile}
        echo -e \"==== (CopyFile) Start Time: \`date\`\"
        cp \${rootPath} \${tagetPath}
        echo -e \"==== (CopyFile) End Time: \`date\`\\n\\n\"
        if [ -f \${tagetPath} ]; then
            echo "Success."
            rm \${rootPath}
            echo \"\${logID}\" >> \${procFile}
        else
            echo "Failure."
        fi
    fi
    
    echo -e \"==== (Job) End Time: \`date\`\"
done

echo -e \"**************************\"
echo -e \"****  FINISH RUNNING  ****\"
echo -e \"**************************\"
" > $job_script
chmod 755 $job_script

unique_job_script=${jobdir}/job__${PROJPATH}__${PROJVERSION}__${proj_title}.sh
cp $job_script $unique_job_script


args_file=${jobdir}/arguments
if [ -f ${args_file} ]; then
    /bin/rm ${args_file}
fi
touch ${args_file}


submit_sub=${jobdir}/submit.sub
echo "executable  = ${unique_job_script}
output      = ${jobdir}/out/CLS\$(ClusterId).PROC\$(ProcId).out
error       = ${jobdir}/err/CLS\$(ClusterId).PROC\$(ProcId).err
log         = ${jobdir}/log/CLS\$(ClusterId).PROC\$(ProcId).log
+JobFlavour = \"${queue}\"
transfer_output_files = \"\"
queue arguments from ${args_file}
" > $submit_sub


submit_script=${jobdir}/submit.sh
echo "#!/bin/bash
echo \"**********************************************************************\"
echo \"****************************** SUBMIT ********************************\"
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
mkdir -p \${dataDir}
if [ ! -d \${dataDir} ]; then
    echo -e \"taget data/ is not exist.\"
    exit
else
    if [ \"\${runmode}\" == \"RUN\" ]; then
        rm -rf \${dataDir}
        mkdir -p \${dataDir}
        if [ ! -d \${dataDir} ]; then
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

procFile=\${jobDir}/proc
if [ ! -f \${procFile} ]; then
    echo -e \"proc is not exist.\"
    exit
fi

if [ \${runmode} == "RUN" ]; then
    /bin/rm \${procFile}
    touch \${procFile}
fi

argsFile=${args_file}
if [ ! -f \${argsFile} ]; then
    echo -e \"argsFile is not exist.\"
    exit
fi

submitSub=${submit_sub}
if [ ! -f \${submitSub} ]; then
    echo -e \"submitSub is not exist.\"
    exit
fi

jobScript=${unique_job_script}
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
        logID=\$(printf \"JOB%07i_EXE%07i\" \${jobID} \${exeID})
        if [ \${runmode} == "RUN" ]; then
            runIt=1
        else
            if [ \"\`grep -n \${logID} \${procFile}\`\" == \"\" ]; then 
                runIt=1
            fi
        fi
    done
    
    if (( \${runIt} == 0 )); then
        continue
    fi

    vargs=\$(printf \"%7i %7i %7i\" \${jobID} \${exeSatID} \${exeEndID})
    echo \"\${vargs}\" >> \${argsFile}

    jobLogID=\$(printf \"ARGS:  JOB %07i    EXE %07i %07i\" \${jobID} \${exeSatID} \${exeEndID})
    echo \"\${jobLogID}\"
done
echo \"***************************** CONDOR *********************************\"
echo \"\"
echo \"COMD:  condor_submit \${submitSub}\"
condor_submit \${submitSub}
echo \"\"
echo \"**********************************************************************\"
echo \"************************** SUBMIT FINISH *****************************\"
echo \"**********************************************************************\"
" >> $submit_script
chmod 755 $submit_script

cfg_clear
