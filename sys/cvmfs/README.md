Hi, Hsin-Yi,
No problem, just add following lines:

CVMFS_REPOSITORIES=ams.cern.ch

CVMFS_HTTP_PROXY="http://csquid01.grid.sinica.edu.tw:3128|http://squid03.grid.sinica.edu.tw:3128|http://squid04.grid.sinica.edu.tw:3128"


Also, since NCU is in Taiwan, it'd better add our repository in domain configuration, so, please ensure "http://cvmfs02.grid.sinica.edu.tw/cvmfs/@fqrn@" is at /etc/cvmfs/domain.d/cern.ch.conf, e.g.:

CVMFS_SERVER_URL="http://cvmfs-stratum-one.cern.ch/cvmfs/@fqrn@;http://cernvmfs.gridpp.rl.ac.uk/cvmfs/@fqrn@;http://cvmfs.racf.bnl.gov/cvmfs/@fqrn@;http://cvmfs.fnal.gov/cvmfs/@fqrn@;http://cvmfs02.grid.sinica.edu.tw/cvmfs/@fqrn@"
CVMFS_KEYS_DIR=/etc/cvmfs/keys/cern.ch
CVMFS_USE_GEOAPI=yes


Let me know if having any questions.

Best regards,
Felix Lee ~
