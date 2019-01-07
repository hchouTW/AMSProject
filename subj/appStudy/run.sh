#!/bin/bash

Version=vdev
Version=test
RunFile=${AMSCore}/subj/appStudy/${Version}/fill

CurDir=${PWD}

DataType=ISS
Stream=${CurDir}/lst/flist.ncu.iss.pass7_19Jan05

OutputDir=${CurDir}/dat

GroupSize=30
Nseq=150

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log
ls -alh /tmp/* | awk '{print $8}' | grep track_fill | xargs /bin/rm" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
