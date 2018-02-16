#!/bin/bash
RunFile=${AMSCore}/subj/appStudy/vdev/app_fill

CurDir=${PWD}

DataType=ISS
Stream=${CurDir}/lst/flist.ncu.iss.pass6_18Feb13

OutputDir=${CurDir}/dat

GroupSize=50
Nseq=60

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}" | qsub -q ams -N JOB${id} -j oe
#    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
