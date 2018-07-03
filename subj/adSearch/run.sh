#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=vdev
#Version=test
RunFile=${AMSCore}/subj/adSearch/${Version}/fill

CurDir=${PWD}

#DataType=ISS
#Stream=${CurDir}/lst/flist.ncu.iss.pass7_18Jun18
DataType=MC
#Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18Jun23
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18Jul03

OutputDir=${CurDir}/dat

#Fit
GroupSize=10
Nseq=60
#Fit
#GroupSize=16
#Nseq=200

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
