#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/vdev/prop_ams02_fill
RunFile=${AMSCore}/subj/trsys/vdev/fit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_smc_fill
DataType=MC
Stream=lst/flist.ncu.mc.pr054000_17Dec12
OutputDir=dat
GroupId=0
#GroupSize=3
GroupSize=35

#${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}

for id in `seq 0 20`
do
    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done
