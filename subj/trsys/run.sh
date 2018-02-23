#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_ams02_fill
RunFile=${AMSCore}/subj/trsys/vdev/fit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_smc_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1119_18Feb05cv
#Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1119_18Jan29
#Stream=${CurDir}/lst/flist.ncu.mc.PR1800_B1042_18Jan29
#Stream=${CurDir}/lst/flist.ncu.mc.PR10004000_B1103_18Jan29

OutputDir=${CurDir}/dat

#Hit
#GroupSize=18
#Nseq=40

#Fit
#GroupSize=9
#Nseq=80
GroupSize=4
Nseq=180

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}" | qsub -q ams -N JOB${id} -j oe
#    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
