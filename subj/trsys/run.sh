#!/bin/bash
Version=vdev
Version=test
#RunFile=${AMSCore}/subj/trsys/${Version}/hit_fill
#RunFile=${AMSCore}/subj/trsys/${Version}/prop_fill
#RunFile=${AMSCore}/subj/trsys/${Version}/bta_fill
RunFile=${AMSCore}/subj/trsys/${Version}/track_fill
#RunFile=${AMSCore}/subj/trsys/${Version}/mu_fill

CurDir=${PWD}

#DataType=ISS
#Stream=${CurDir}/lst/flist.ncu.iss.pass7_19Jan08
DataType=MC
Stream=${CurDir}/lst/flist.ncu.mc.HE4_250T_B1200_19Jan09
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep21
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_24000_B1200_18Sep21
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Oct17
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Dec23
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_24000_B1200_18Oct17
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_24000_B1200_18Dec23
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_216000_B1200_18Dec23
#Stream=${CurDir}/lst/flist.ncu.mc.EL_025500_B1200_18Oct03

OutputDir=${CurDir}/dat

#Hit
#GroupSize=20
#Nseq=230
#Fit
GroupSize=15
Nseq=300

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log
ls -alh /tmp/* | awk '{print $8}' | grep fill | xargs /bin/rm" | qsub -q ams -N JOB${id} -j oe
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done

#RunFile=${AMSCore}/subj/trsys/${Version}/hitq_fit
#for id in `seq 1 100`
#do
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#${RunFile} ${id}" | qsub -q ams -N JOB${id} -j oe
#done

#/opt/glibc-2.17/lib/ld-2.17.so --library-path /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6/lib64:/opt/glibc-2.17/lib:/usr/lib64:/lib64:$LD_LIBRARY_PATH /ams_home/hchou/AMSCore/subj/trsys/vdev/track_fill MC lst/flist.cern.mc.PR054000_B1200_18Jul04 0 1
