#!/bin/bash
Proj=${AMSCore}/subj/ccStudy/vdev

RunFile=${Proj}/cc_fill

DataType=MC
Stream=${Proj}/lst/flist.ncu.mc.PR054000_B1119_18Feb05cv2

OutputDir=${PWD}/dat

#Hit
#GroupSize=18
#Nseq=40

#Fit
#GroupSize=9
#Nseq=80
GroupSize=12
Nseq=65
#Nseq=13

for id in `seq 0 ${Nseq}`
do
    echo "%!/bin/bash
source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}" | qsub -q ams -N JOB${id}
#    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
