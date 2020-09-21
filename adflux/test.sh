#/bin/bash
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_acc
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_rlt
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_app
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_tme
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/phys_tme

#exe/mdst_sel -type=ISS -inpath=lst/flist.cern.iss.B1130.pass7 -gindex=633 -gsize=50
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.021000.B1220 -gindex=33 -gsize=5
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.apfluxr.pl1.l1.021000.B1220 -gindex=32 -gsize=10
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.flux.504000.B1220 -gindex=0  -gsize=1
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.l19.10016000.B1220 -gindex=0  -gsize=100
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.l1o9.flux.2016000.B1220 -gindex=0  -gsize=10
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=1  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=2  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=3  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=4  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=5  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=6  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=7  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=8  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=9  -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=10 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=11 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=12 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=13 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=14 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=15 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=16 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=17 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=18 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=19 -gsize=1 &
#exe/mdst_sel -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=20 -gsize=1 &

exe/mdst_ana -type=ISS -inpath=lst/flist.cern.proj.iss.pass7 -gindex=3 -gsize=1
#exe/mdst_ana -type=MC -inpath=lst/flist.cern.proj.mc.prL -gindex=3 -gsize=1

#for idx in {0..35}
#do
#    exe/lvtme ${idx} &> /dev/null &
#done
