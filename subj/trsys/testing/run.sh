#!/bin/bash
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop_ams02_fill
#RunFile=${AMSCore}/subj/trsys/testing/vdev/fit_ams02_fill
DataType=MC
Stream=lst/flist.ncu.mc.pr0510_17Oct30
#Stream=lst/flist.ncu.mc.pr1800_17Oct30
#Stream=lst/flist.ncu.mc.pr05800_17Oct30
#Stream=lst/flist.mc.pr2016000
OutputDir=dat
GroupId=0
GroupSize=5

#${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}

for id in `seq 0 20`
do
    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
