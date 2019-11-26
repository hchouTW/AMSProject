#/bin/bash
#mkdir -p out/doc

#exe/fit_mscat    19Jul08/test16_pr
#exe/fit_eloss    19Jul08/test16_pr
#exe/fit_geom     19Jul08/test16_pr
#exe/fit_vel_tf   19Jul08/test16_pr
#exe/fit_vel_rh   19Jul08/test16_pr
#exe/fit_phys_tf  19Jul08/test16_pr
#exe/fit_phys_rh  19Jul08/test16_pr
#exe/fit_mutr_tf  19Jul08/test16_pr
#exe/fit_mutr_rh  19Jul08/test16_pr

#exe/fit_mscat    19Jul08/test16_he
#exe/fit_eloss    19Jul08/test16_he
#exe/fit_geom     19Jul08/test16_he
#exe/fit_vel_tf   19Jul08/test16_he
#exe/fit_vel_rh   19Jul08/test16_he
#exe/fit_phys_tf  19Jul08/test16_he
#exe/fit_phys_rh  19Jul08/test16_he
#exe/fit_mutr_tf  19Jul08/test16_he
#exe/fit_mutr_rh  19Jul08/test16_he

exe/fit_merge out/doc_pr5 out/doc_he5

#exe/mdst -type=MC -inpath=lst/flist.cern.mc.he4.pl1.l19.216000.B1200 -gindex=60 -gsize=1
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=60 -gsize=1
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=0  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=1  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=2  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=3  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=4  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=5  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=6  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=7  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=8  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=9  -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=10 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=11 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=12 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=13 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=14 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=15 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=16 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=17 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=18 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=19 -gsize=1 &
#exe/mdst -type=MC -inpath=lst/flist.cern.mc.pr.pl1.l1.054000.B1200 -gindex=20 -gsize=1 &
