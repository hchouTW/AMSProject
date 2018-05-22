#!/bin/bash

for id in `seq 0 22`
do
    ii=`printf "%05d" $id`
    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May19/PR054000.B1200 /ams02/user/hchou/PR054000.B1200_18May19 YiNtuple_MC.${ii} &
done
