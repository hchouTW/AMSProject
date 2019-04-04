#!/bin/bash
Version=vdev
Version=19Apr04

ClassDef=${AMSCore}/prod/${Version}/lib
LD_LIBRARY_PATH=${ClassDef}:${LD_LIBRARY_PATH}
RunFile=${AMSCore}/prod/${Version}/YiProdNtuple

CurDir=${PWD}

DataType=ISS
Stream=${CurDir}/lst/flist.cern.iss.B1130.pass7
#Stream=${CurDir}/lst/flist.cern.iss.B1130.pass7.test
#Stream=${CurDir}/lst/flist.ncu.iss.B1130.pass7
#Stream=${CurDir}/lst/flist.loc.iss.B1130.pass7
#DataType=MC
#Stream=${CurDir}/lst/flist.cern.mc.el.pl1.025200.B1200
#Stream=${CurDir}/lst/flist.cern.mc.el.pl1.2004000.B1200
#Stream=${CurDir}/lst/flist.cern.mc.pr.pl1.l1.054000.B1200
#Stream=${CurDir}/lst/flist.cern.mc.he4.pl1.l1.24000.B1200
#Stream=${CurDir}/lst/flist.cern.mc.he4.pl1.l19.216000.B1200
#Stream=${CurDir}/lst/flist.cern.mc.he4.pl1.l1.150T.B1200
#DataType=BT
#Stream=${CurDir}/lst/flist.cern.bt.pr.400.B1082
#Stream=${CurDir}/lst/flist.ncu.bt.pr.400.B1082

GroupId=5
#GroupId=19912
#GroupId=5005
GroupSize=1
OutputDir=.

if [[ -f ${RunFile} && -f ${Stream} ]]; then
  ${RunFile} ${DataType} ${Stream} ${GroupId} ${GroupSize} ${OutputDir}
else
  echo -e "Please Check File :"
  if [ ! -f ${RunFile} ]; then echo ${RunFile}; fi
  if [ ! -f ${Stream} ]; then echo ${Stream}; fi
fi


#Nseq=81
#GroupSize=1
#OutputDir=${CurDir}/dat2
#for id in `seq 0 ${Nseq}`
#do
#    echo "%!/bin/bash
#source /ams_home/hchou/AMSProject/sw/ROOT/setup_amsenv_root5gcc.sh
#export LD_LIBRARY_PATH=${ClassDef}:\${LD_LIBRARY_PATH}
#echo \"${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir}\" > ${CurDir}/log/JOB${id}.log
#${RunFile} ${DataType} ${Stream} ${id} ${GroupSize} ${OutputDir} >> ${CurDir}/log/JOB${id}.log" | qsub -q ams -N JOB${id} -j oe
#done
