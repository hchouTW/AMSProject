#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=vdev
RunFile=${AMSCore}/subj/adSearch/${Version}/fill

CurDir=${PWD}

#DataType=ISS
#Stream=${CurDir}/lst/flist.ncu.iss.pass7.mfixed_B1130_18May15
DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18May17

OutputDir=${CurDir}/dat

GroupSize=4
Nseq=50
#Fit
#GroupSize=3
#Nseq=450
#GroupSize=10
#Nseq=150

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
done
