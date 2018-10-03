#!/bin/bash
Version=vdev
Version=18Oct03

ClassDef=${AMSCore}/prod/${Version}/lib
LD_LIBRARY_PATH=${ClassDef}:${LD_LIBRARY_PATH}
RunFile=${AMSCore}/prod/${Version}/YiProdNtuple

#DataType=ISS
#Stream=lst/flist.cern.iss.B1130.pass7
#Stream=lst/flist.ncu.iss.B1130.pass7
DataType=MC
Stream=lst/flist.cern.mc.pr.pl1.l1.054000.B1200
#Stream=lst/flist.cern.mc.he4.pl1.l1.24000.B1200
#DataType=BT
#Stream=lst/flist.cern.bt.pr.400.B1082

GroupId=10
GroupSize=1
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
