#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill

Version=vdev
Version=test
RunFile=${AMSCore}/subj/trsys/${Version}/hit_fill
#RunFile=${AMSCore}/subj/trsys/${Version}/prop_fill
#RunFile=${AMSCore}/subj/trsys/${Version}/track_fill

CurDir=${PWD}

DataType=MC
#Stream=/ams_home/hchou/tmp/flist.pr
#Stream=${CurDir}/lst/flist.ncu.mc.PR054000_B1200_18Jul04
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep16
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_24000_B1200_18Sep17
#Stream=${CurDir}/lst/flist.ncu.mc.C12_612000_B1200_18Sep18
#Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep20
#Stream=${CurDir}/lst/flist.ncu.mc.HE4_24000_B1200_18Sep21
Stream=${CurDir}/lst/flist.ncu.mc.PR_054000_B1200_18Sep21

OutputDir=${CurDir}/dat

#Hit
GroupSize=75
Nseq=39
#Fit
#GroupSize=6
#Nseq=500

for id in `seq 0 ${Nseq}`
do
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} > ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done

#RunFile=${AMSCore}/subj/trsys/${Version}/hit_fit
#for id in `seq 15 150`
#do
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#${RunFile} ${id}" | qsub -q ams -N JOB${id} -j oe
#done


#/opt/glibc-2.17/lib/ld-2.17.so --library-path /cvmfs/sft.cern.ch/lcg/external/gcc/6.2.0/x86_64-slc6/lib64:/opt/glibc-2.17/lib:/usr/lib64:/lib64:$LD_LIBRARY_PATH /ams_home/hchou/AMSCore/subj/trsys/vdev/track_fill MC lst/flist.cern.mc.PR054000_B1200_18Jul04 0 1
