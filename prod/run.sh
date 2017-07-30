#!bin/bash
YiBinDir=/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/core/production
YiDataSet=/afs/cern.ch/user/h/hchou/public/BSUB/submit/user/hchou/dataset/production

RunFile=${YiBinDir}/vdev/YiProdNtuple
#RunFile=${YiBinDir}/V.2017Apr03/YiProdNtuple
#DataType=ISS
DataType=MC
#Stream=${YiDataSet}/iss.B950.pass6
#Stream=${YiDataSet}/mc.ap.pl1.0_510.B1042
#Stream=${YiDataSet}/mc.ap.pl1.1800.B1042
#Stream=${YiDataSet}/mc.pr.pl1.0_510.B1042
#Stream=${YiDataSet}/mc.pr.pl1.1800.B1042
#Stream=${YiDataSet}/mc.pr.pl1.1200.c090.B1030
#Stream=${YiDataSet}/mc.pr.pl1.1200.c110.B1030
Stream=${YiDataSet}/mc.pr.pl1.flux.l1a9.2016000.B1042
#Stream=${YiDataSet}/mc.pr.pl1.flux.l1o9.2016000.B1042
#Stream=${YiDataSet}/mc.el.pl1.0_25200.B1086
#Stream=${YiDataSet}/mc.el.pl1.2002000.B1086
#Stream=${YiDataSet}/mc.ph.pl1.0_050_25.B1069
#Stream=${YiDataSet}/mc.ph.pl1.0_2510.B1069
#Stream=${YiDataSet}/mc.ph.pl1.101000.B1069
GroupId=250
GroupSize=2
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi
