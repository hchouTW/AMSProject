#!/bin/bash
#RunFile=${AMSCore}/subj/trsys/testing/vdev/prop
#DataType=MC
#Stream=lst/flist.proton6
##Stream=lst/flist.proton8
#OutputDir=.
#GroupId=1
#GroupSize=0

#RunFile=${AMSCore}/subj/trsys/testing/vdev/fit
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop_ams02
DataType=MC
Stream=lst/flist.ncu.mc.pr1800
#Stream=lst/flist.ncu.mc.pr0510
#Stream=lst/flist.mc.pr2016000
OutputDir=.
GroupId=0
GroupSize=0

${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
