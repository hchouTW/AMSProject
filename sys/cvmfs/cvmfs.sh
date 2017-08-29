yum localinstall http://cvmrepo.web.cern.ch/cvmrepo/yum/cvmfs/EL/6/x86_64/cvmfs-release-2-4.el6.noarch.rpm -y
yum install cvmfs cvmfs-init-scripts cvmfs-auto-setup -y
echo CVMFS_HTTP_PROXY=http://ca-proxy.cern.ch:3128 >> /etc/cvmfs/default.local
ls /cvmfs/ams.cern.ch/Offline
