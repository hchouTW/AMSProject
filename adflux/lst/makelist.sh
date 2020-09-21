#!/bin/bash
alias readflist=${AMSProjJobs}/CERN/readflist.sh
subv=33
readflist flist.cern.proj.iss.pass7      /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/iss${subv}
readflist flist.cern.proj.mc.pr_l1o9flux /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcpr_l1o9flux${subv}
readflist flist.cern.proj.mc.ap          /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcap${subv}
readflist flist.cern.proj.mc.pr          /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcpr${subv}
readflist flist.cern.proj.mc.prL         /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcprL${subv}
readflist flist.cern.proj.mc.ad          /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcad${subv}
readflist flist.cern.proj.mc.d           /eos/ams/user/h/hchou/AMSData/proj/adflux/20Jan15/mcd${subv}
