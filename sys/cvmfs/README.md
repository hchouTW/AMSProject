## Install CVMFS
In SLC6 system
```
yum install cvmfs.x86_64
```
In other system
```
yum localinstall http://cvmrepo.web.cern.ch/cvmrepo/yum/cvmfs/EL/6/x86_64/cvmfs-release-2-4.el6.noarch.rpm -y
yum install cvmfs cvmfs-init-scripts cvmfs-auto-setup -y
```
## CVMFS Setting

In /etc/cvmfs/default.local file
```
CVMFS_REPOSITORIES=ams.cern.ch
CVMFS_DEFAULT_DOMAIN=cern.ch
CERNVM_SERVER_URL="http://cvmfsrep.grid.sinica.edu.tw/opt/@org@;http://cvmfs-stratum-one.cern.ch/opt/@org@;http://cernvmfs.gridpp.rl.ac.uk/opt/@org@;http://cvmfs.racf.bnl.gov/opt/@org@"
CVMFS_SERVER_URL="http://cvmfsrep.grid.sinica.edu.tw/opt/@org@;http://cvmfs-stratum-one.cern.ch/opt/@org@;http://cernvmfs.gridpp.rl.ac.uk/opt/@org@;http://cvmfs.racf.bnl.gov/opt/@org@"
CVMFS_HTTP_PROXY="http://csquid01.grid.sinica.edu.tw:3128|http://squid03.grid.sinica.edu.tw:3128|http://squid04.grid.sinica.edu.tw:3128"
CVMFS_QUOTA_LIMIT=11000
CVMFS_MEMCACHE_SIZE=48
```

In /etc/cvmfs/domain.d/cern.ch.conf file
```
CVMFS_SERVER_URL="http://cvmfs-stratum-one.cern.ch/cvmfs/@fqrn@;http://cernvmfs.gridpp.rl.ac.uk/cvmfs/@fqrn@;http://cvmfs.racf.bnl.gov/cvmfs/@fqrn@;http://cvmfs.fnal.gov/cvmfs/@fqrn@;http://cvmfs02.grid.sinica.edu.tw/cvmfs/@fqrn@"
CVMFS_KEYS_DIR=/etc/cvmfs/keys/cern.ch
CVMFS_USE_GEOAPI=yes
```

In /etc/fuse.conf file
```
user_allow_other # added by CernVM-FS
```

In /etc/auto.master file
```
+auto.master
/cvmfs program:/etc/auto.cvmfs
```

Restart autofs Service
```
service autofs restart
```
