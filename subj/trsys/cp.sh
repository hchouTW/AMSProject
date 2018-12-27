#!/bin/bash

#for id in `seq 0 33`
for id in `seq 0 41`
do
    ii=`printf "%05d" $id`
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct17/MC_EL_025500.B1200 /ams02/user/hchou/MC_EL_025500.B1200_18Oct17 YiNtuple_MC.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Dec23/MC_PR_054000.B1200 /ams02/user/hchou/MC_PR_054000.B1200_18Dec23 YiNtuple_MC.${ii} &
    sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Dec23/MC_HE4_24000.B1200 /ams02/user/hchou/MC_HE4_24000.B1200_18Dec23 YiNtuple_MC.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct11/BT_PR_400.B1082 /ams02/user/hchou/BT_PR_400.B1082_18Oct11 YiNtuple_BT.${ii} &
    #sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18Oct11/BT_EL_180.B1082 /ams02/user/hchou/BT_EL_180.B1082_18Oct11 YiNtuple_BT.${ii} &
done

#sh /ams_home/hchou/AMSProject/jobs/NCU/xrdcp_cern2ncu.sh /eos/ams/user/h/hchou/AMSData/prod/18May23/CERN_TEST /ams02/user/hchou/iss.pass7.B1130_18May23 YiNtuple_ISS &
