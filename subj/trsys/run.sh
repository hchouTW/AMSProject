#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/landau
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1119_18Feb27
#Stream=${CurDir}/lst/flist.ncu.mc.PR1800_B1042_18Jan29
#Stream=${CurDir}/lst/flist.ncu.mc.PR10004000_B1103_18Jan29

OutputDir=${CurDir}/dat

#Hit
#GroupSize=15
#Nseq=100
#Fit
GroupSize=6
Nseq=250

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#${RunFile} ${id} > ${CurDir}/log/JOB${id}" | qsub -q ams -N JOB${id} -j oe
done
