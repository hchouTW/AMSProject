# Install castor-lib.x86_64 (AMS software requirement)

## In SLC6 system
```bash
yum install castor-lib.x86_64
```

## In CentOS6 system
In /etc/yum.repos.d/slc6-extras.repo file
```
[slc6-extras]
name=Scientific Linux CERN 6 (SLC6) add-on packages, no formal support
baseurl=http://linuxsoft.cern.ch/cern/slc6X/$basearch/yum/extras/
gpgkey=file:///etc/pki/rpm-gpg/RPM-GPG-KEY-cern
gpgcheck=0
enabled=1
protect=1
priority=5
```
check
```bash
yum search castor-lib
```
install
```bash
yum install castor-lib.x86_64
```
