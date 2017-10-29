#!/bin/bash
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop_ams02_fill
DataType=MC
#Stream=lst/flist.ncu.mc.pr05800
Stream=lst/flist.mc.pr2016000
OutputDir=dat
GroupId=0
GroupSize=7

#${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}

for id in `seq 0 9`
do
    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
