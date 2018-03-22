#!/bin/bash

for id in `seq 0 24`
do
    ii=`printf "%05d" $id`
    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Mar22/PR054000.B1119 /ams02/user/hchou/PR054000.B1119_18Mar22 YiNtuple_MC.${ii} &
done
