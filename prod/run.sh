#!/bin/bash
Version=vdev

ClassDef=${AMSCore}/prod/${Version}/lib
LD_LIBRARY_PATH=${ClassDef}:${LD_LIBRARY_PATH}
RunFile=${AMSCore}/prod/${Version}/YiProdNtuple

DataType=ISS
Stream=lst/flist.ncu.iss.B950.pass6
#DataType=MC
#Stream=lst/flist.cern.mc.pr.pl1.0510
#Stream=lst/flist.cern.mc.pr.pl1.1800

GroupId=1
GroupSize=1
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
