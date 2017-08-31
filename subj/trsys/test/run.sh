#!bin/bash
RunFile=${AMSCore}/subj/trsys/test/testing
DataType=ISS
Stream=list/flist.iss.B950.pass6
OutputDir=.
GroupId=5
GroupSize=6
#${RunFile} -t ${DataType} -i ${Stream} -o ${OutputDir} -g ${GroupId} -s ${GroupSize}
${RunFile} -t ${DataType} -i ${Stream} -o ${OutputDir} -g ${GroupId} -s ${GroupSize} 








#
#if [[ -f ${RunFile} && -f ${Stream} ]]; then
#  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
#else
#  echo -e "Please Check File :"
#  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
#  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
#fi
