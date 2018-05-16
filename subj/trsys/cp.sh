#!/bin/bash

#for id in `seq 0 16`
#do
#    ii=`printf "%05d" $id`
#    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May15/PR054000.B1200 /ams02/user/hchou/PR054000.B1200_18May15 YiNtuple_MC.${ii} &
#done

sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May17/PR054000.B1200 /ams02/user/hchou/PR054000.B1200_18May17 YiNtuple_MC* &
