#!/bin/bash
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop
DataType=ISS
#Stream=lst/flist.proton6
Stream=lst/flist.proton8
OutputDir=.
GroupId=1
GroupSize=0

${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
