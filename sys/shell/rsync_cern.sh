#!/usr/bin/expect -f
set password [lindex $argv 0]
set inpath [lindex $argv 1]
set outpath [lindex $argv 2]
spawn rsync --ignore-existing --progress --human-readable hchou@lxplus.cern.ch:$inpath $outpath
expect {
    "*Password:" { send "$password\r"}
}
interact
