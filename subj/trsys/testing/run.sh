#!/bin/bash
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop
DataType=ISS
Stream=lst/flist.proton2
OutputDir=.
GroupId=0
GroupSize=0

${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
