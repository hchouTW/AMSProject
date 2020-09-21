#!/bin/bash
#sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.iss.B1130.pass7
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.pr.l1o9.flux.2016000.B1220
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.ap.pl1.021000.B1220 
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.pr.pl1.021000.B1220 
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.pr.pl1.0550.B1220 
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.ad.pl1.021000.B1128 
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.mc.d.pl1.021000.B1128 
