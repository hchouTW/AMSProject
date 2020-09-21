#!/bin/bash
sh ${AMSProjJobs}/CERN/submit_condor.sh RUN jobconf.cern.iss.B1130.pass7
sh ${AMSProjJobs}/CERN/submit_condor.sh RUN jobconf.cern.mc.ap.pl1.l1.021000.B1220 
sh ${AMSProjJobs}/CERN/submit_condor.sh RUN jobconf.cern.mc.d.pl1.l1.021000.B1128 
sh ${AMSProjJobs}/CERN/submit_condor.sh RUN jobconf.cern.mc.pr.pl1.l1.021000.B1220 
sh ${AMSProjJobs}/CERN/submit_condor.sh RUN jobconf.cern.mc.pr.0550.B1220 
