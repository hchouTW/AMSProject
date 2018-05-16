#!/bin/bash

#for id in `seq 0 16`
#do
#    ii=`printf "%05d" $id`
#    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May15/PR054000.B1200 /ams02/user/hchou/PR054000.B1200_18May15 YiNtuple_MC.${ii} &
#done

sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May15/CERN_TEST /ams02/user/hchou/iss_pass7_18May15 YiNtuple_ISS* &
