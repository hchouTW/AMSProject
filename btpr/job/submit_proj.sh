#!/bin/bash
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.iss.pass7
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.pr.l1o9flux
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.pr
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.prL
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.ap
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.d
sh ${AMSProjJobs}/CERN/submit_condor.sh RERUN jobconf.cern.proj.mc.ad