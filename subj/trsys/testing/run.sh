#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/testing/vdev/prop
#DataType=MC
#Stream=lst/flist.proton6
#Stream=lst/flist.proton8
#OutputDir=.
#GroupId=1
#GroupSize=0

RunFile=${AMSCore}/subj/trsys/testing/vdev/fit
DataType=MC
Stream=lst/flist.mc.pr2016000
OutputDir=.
GroupId=1
GroupSize=15

${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
