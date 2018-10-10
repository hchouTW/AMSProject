#!/bin/bash

for id in `seq 0 33`
#for id in `seq 0 41`
do
    ii=`printf "%05d" $id`
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct03/EL_025500.B1200 /ams02/user/hchou/EL_025500.B1200_18Oct03 YiNtuple_MC.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct03/PR_054000.B1200 /ams02/user/hchou/PR_054000.B1200_18Oct03 YiNtuple_MC.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Sep21/HE4_24000.B1200 /ams02/user/hchou/HE4_24000.B1200_18Sep21 YiNtuple_MC.${ii} &
    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct10/PR_400.B1082 /ams02/user/hchou/BT_PR_400.B1082_18Oct10 YiNtuple_BT.${ii} &
done

#sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May23/CERN_TEST /ams02/user/hchou/iss.pass7.B1130_18May23 YiNtuple_ISS &
