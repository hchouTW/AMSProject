#/bin/bash
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_acc
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_rlt
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_app
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/adflux_tme
#/afs/cern.ch/user/h/hchou/AMSProject/adflux/exe/phys_tme

#exe/mdst_sel -type=BT -inpath=lst/flist.cern.bt.pr.400.B1082 -gindex=3 -gsize=1

exe/mdst_ana -type=BT -inpath=lst/flist.cern.proj.bt.pr -gindex=0 -gsize=1

#for idx in {0..35}
#do
#    exe/lvtme ${idx} &> /dev/null &
#done
