#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_ams02_fill
RunFile=${AMSCore}/subj/trsys/vdev/fit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_smc_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.pr054000_17Dec23
#Stream=${CurDir}/lst/flist.ncu.mc.el025200_17Dec23
#Stream=${CurDir}/lst/flist.ncu.mc.el2004000_17Dec23

OutputDir=${CurDir}/dat

#Hit
#GroupSize=18
#Nseq=40

#Fit
GroupSize=4
Nseq=180
#GroupSize=8
#Nseq=40

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}" | qsub -q ams -N JOB${id} -j oe
#    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
