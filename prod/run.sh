#!/bin/bash
Version=17Nov24

ClassDef=${AMSCore}/prod/${Version}/lib
LD_LIBRARY_PATH=${ClassDef}:${LD_LIBRARY_PATH}
RunFile=${AMSCore}/prod/${Version}/YiProdNtuple

#DataType=ISS
#Stream=lst/flist.ncu.iss.B950.pass6
DataType=MC
#Stream=lst/flist.cern.mc.pr.pl1.0510
#Stream=lst/flist.cern.mc.pr.pl1.1800
#Stream=lst/flist.cern.mc.pr.pl1.flux.l1a9.2016000
#Stream=lst/flist.cern.mc.pr.pl1.l1.054000.B1119
#Stream=lst/flist.cern.mc.pr.pl1.05100.B1119
#Stream=lst/flist.cern.mc.el.pl1.0_25200.B1118
Stream=lst/flist.cern.mc.el.pl1.2004000.B1118

GroupId=2
GroupSize=50
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
