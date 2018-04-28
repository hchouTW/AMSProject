#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/time_fill
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
#RunFile=${AMSCore}/subj/trsys/vdev/track_fill
RunFile=${AMSCore}/subj/trsys/test8/track_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.cern.mc.PR054000_B1119_18Mar23
#Stream=${CurDir}/lst/flist.cern.mc.EL025200_B1119_18Mar23

OutputDir=/afs/cern.ch/work/h/hchou/AMSData/test41
mkdir -p $OutputDir
mkdir -p $OutputDir/log

#Hit
#GroupSize=30
#Nseq=50
#FitGroupSize=8
#GroupSize=5
#Nseq=250
GroupSize=1
Nseq=1400

#for id in `seq 1 100`
for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /afs/cern.ch/user/h/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
source /eos/ams/user/h/hchou/ExternalLibs/BIN/CERN_TRACKSys.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}" | bsub -q ams1nd -J JOB${id} -oo ${OutputDir}/log/JOB${id}.log
#${RunFile} ${id}" | bsub -q ams1nd -J JOB${id} -oo ${OutputDir}/log/JOB${id}.log
done
