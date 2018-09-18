#!/bin/bash

for id in `seq 0 39`
#for id in `seq 350 370`
do
    ii=`printf "%05d" $id`
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Sep16/PR_054000.B1200 /ams02/user/hchou/PR_054000.B1200_18Sep16 YiNtuple_MC.${ii} &
    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Sep17/HE4_24000.B1200 /ams02/user/hchou/HE4_24000.B1200_18Sep17 YiNtuple_MC.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Jun18/CERN_TEST /ams02/user/hchou/iss.pass7_18Jun18 YiNtuple_ISS.${ii} &
done

#sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May23/CERN_TEST /ams02/user/hchou/iss.pass7.B1130_18May23 YiNtuple_ISS &
