#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=test
#RunFile=${AMSCore}/subj/trsys/${Version}/hit_fill
RunFile=${AMSCore}/subj/trsys/${Version}/track_fill

CurDir=${PWD}

DataType=MC
#Stream=/ams_home/hchou/tmp/flist.pr
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18May19

OutputDir=${CurDir}/dat

#Hit
#GroupSize=20
#Nseq=70
#Fit
#GroupSize=3
#Nseq=450
GroupSize=4
Nseq=500

#for id in `seq 1 100`
for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${id} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
done
