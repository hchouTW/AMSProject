#!bin/bash
RunFile=${AMSCore}/subj/template/YiSlcNtuple
DataType=ISS
Stream=lst/flist.iss.B950.pass6
GroupId=0
GroupSize=1
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
