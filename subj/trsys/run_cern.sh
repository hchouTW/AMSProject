#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fit
#RunFile=${AMSCore}/subj/trsys/vdev/hit_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_fill
RunFile=${AMSCore}/subj/trsys/vdev/track_fill

CurDir=${PWD}

DataType=MC
Stream=${CurDir}/lst/flist.cern.mc.PR054000_B1119_18Mar23

OutputDir=/afs/cern.ch/work/h/hchou/AMSData/test11
mkdir -p $OutputDir
mkdir -p $OutputDir/log

#Hit
#GroupSize=30
#Nseq=65
#Fit
GroupSize=10
Nseq=200

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /afs/cern.ch/user/h/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
source /eos/ams/user/h/hchou/ExternalLibs/BIN/CERN_TRACKSys.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}" | bsub -q ams1nd -J JOB${id} -oo ${OutputDir}/log/JOB${id}.log
done
