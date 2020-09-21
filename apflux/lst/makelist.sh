#!/bin/bash
alias readflist=${AMSProjJobs}/CERN/readflist.sh
subv=51
readflist flist.cern.proj.iss.pass7      /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/iss${subv}
readflist flist.cern.proj.mc.ap          /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcap${subv}
readflist flist.cern.proj.mc.ap_patch    /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcap${subv}patch
readflist flist.cern.proj.mc.ap_minus    /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcap_minus${subv}
readflist flist.cern.proj.mc.ap_plus     /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcap_plus${subv}
readflist flist.cern.proj.mc.pr          /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcpr${subv}
readflist flist.cern.proj.mc.pr_l1o9flux /eos/ams/user/h/hchou/AMSData/proj/apflux/20Jan15/mcpr_l1o9flux${subv}
