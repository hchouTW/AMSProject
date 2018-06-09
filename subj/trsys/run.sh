#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=vdev
Version=test
#RunFile=${AMSCore}/subj/trsys/${Version}/hit_fill
RunFile=${AMSCore}/subj/trsys/${Version}/track_fill

CurDir=${PWD}

DataType=MC
#Stream=/ams_home/hchou/tmp/flist.pr
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18May27
#Stream=${CurDir}/lst/flist.ncu.mc.EL2004000_B1119_18May27

OutputDir=${CurDir}/dat

#Hit
#GroupSize=50
#Nseq=40
#GroupSize=24
#Nseq=20
#Fit
GroupSize=7
Nseq=300
#GroupSize=4
#Nseq=500

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
