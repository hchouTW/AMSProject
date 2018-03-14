#!/bin/bash
RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1119_18Mar12

OutputDir=${CurDir}/dat

#Hit
GroupSize=30
Nseq=65
#Fit
#GroupSize=8
#Nseq=120

#for id in `seq 0 ${Nseq}`
for id in `seq 1 150`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${id} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
done
