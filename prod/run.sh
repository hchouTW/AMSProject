#!/bin/bash
Version=vdev
#Version=18Sep18

ClassDef=${AMSCore}/prod/${Version}/lib
LD_LIBRARY_PATH=${ClassDef}:${LD_LIBRARY_PATH}
RunFile=${AMSCore}/prod/${Version}/YiProdNtuple

#DataType=ISS
#Stream=lst/flist.cern.iss.B1130.pass7
#Stream=lst/flist.ncu.iss.B1130.pass7
DataType=MC
#Stream=lst/flist.cern.mc.pr.pl1.l1.054000.B1200
#Stream=lst/flist.cern.mc.he4.pl1.l1.24000.B1200
Stream=lst/flist.cern.mc.c12.pl1.l1.612000.B1200

GroupId=1005
GroupSize=1
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
