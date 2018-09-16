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
#Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18Jul04
Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep16

OutputDir=${CurDir}/dat

#Fit
#GroupSize=5
#Nseq=110
#Fit
GroupSize=20
Nseq=22

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
