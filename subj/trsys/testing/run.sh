#!bin/bash
RunFile=${AMSCore}/subj/trsys/testing/prop
DataType=ISS
Stream=lst/flist.iss.B950.pass6
OutputDir=.
GroupId=5
GroupSize=6
#${RunFile} -t ${DataType} -i ${Stream} -o ${OutputDir} -g ${GroupId} -s ${GroupSize}
#${RunFile} -t ${DataType} -i ${Stream} -o ${OutputDir} -g ${GroupId} -s ${GroupSize} 
${RunFile} ISS lst/flist.proton 








#
#if [[ -f ${RunFile} && -f ${Stream} ]]; then
#  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
#else
#  echo -e "Please Check File :"
#  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
#  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
#fi
