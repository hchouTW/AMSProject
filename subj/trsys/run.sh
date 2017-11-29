#!/bin/bash
RunFile=${AMSCore}/subj/trsys/vdev/prop_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/fit_ams02_fill
#RunFile=${AMSCore}/subj/trsys/vdev/prop_smc_fill
DataType=MC
#Stream=lst/flist.ncu.mc.pr05100_17Nov24
Stream=lst/flist.ncu.mc.pr054000_17Nov24
OutputDir=dat
GroupId=0
#GroupSize=5
GroupSize=5

#${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}

for id in `seq 0 20`
#for id in `seq 0 10`
do
    ${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} &> /dev/null &
done

#ljcheck fill
