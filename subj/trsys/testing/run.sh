#!bin/bash
RunFile=${AMSCore}/subj/trsys/testing/vdev/prop
DataType=ISS
Stream=lst/flist.proton
OutputDir=.
GroupId=5
GroupSize=6
${RunFile} ${DataType} ${Stream}








#
#if [[ -f ${RunFile} && -f ${Stream} ]]; then
#  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
#else
#  echo -e "Please Check File :"
#  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
#  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
#fi
