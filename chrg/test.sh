#/bin/bash
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_acc
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_rlt
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_app
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_tme
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/phys_tme

#exe/mdst_sel -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7 -gindex=633 -gsize=10
#exe/mdst_sel -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7.20Jul01 -gindex=3 -gsize=10
#exe/mdst_sel -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7.20Jul08 -gindex=3 -gsize=1
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.B1220 -gindex=1 -gsize=10
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.fe.B1220 -gindex=0 -gsize=1
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.c.B1215 -gindex=0 -gsize=1

#exe/mdst_ana -type=ISS -inpath=lst/flist.cern.proj.iss.pass7 -gindex=0 -gsize=1
#exe/mdst_ana -type=MC -inpath=lst/flist.cern.proj.mc.pr -gindex=3 -gsize=1
#exe/mdst_ana -type=MC -inpath=lst/flist.cern.proj.mc.c -gindex=3 -gsize=1
#exe/mdst_ana -type=MC -inpath=lst/flist.cern.proj.mc.fe -gindex=3 -gsize=1

#for idx in {0..35}
#do
#    exe/lvtme ${idx} &> /dev/null &
#done
