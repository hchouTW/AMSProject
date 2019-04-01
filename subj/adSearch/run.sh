#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=vdev
Version=test
RunFile=${AMSCore}/subj/adSearch/${Version}/fill

CurDir=${PWD}

DataType=ISS
Stream=${CurDir}/lst/flist.ncu.iss.rich
Stream=${CurDir}/lst/flist.ncu.iss.rich5
#Stream=${CurDir}/lst/flist.ncu.iss.pass7_18Jun18
#Stream=${CurDir}/lst/flist.ncu.iss.pass7_18Dec23
#Stream=${CurDir}/lst/flist.ncu.iss.pass7_19Jan05
#DataType=MC
#Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18Jul04
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep16
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Oct17
#Stream=${CurDir}/lst/flist.ncu.mc.EL_025500_B1200_18Oct03
#DataType=BT
#Stream=${CurDir}/lst/flist.ncu.bt.PR_400_B1082_18Oct10

OutputDir=${CurDir}/dat

#Fit
GroupSize=35
Nseq=30
#Fit
#GroupSize=1
#Nseq=52

for id in `seq 0 ${Nseq}`
do
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log
#ls -alh /tmp/* | awk '{print $8}' | grep hchou | xargs /bin/rm" | qsub -q ams -N JOB${id} -j oe
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
