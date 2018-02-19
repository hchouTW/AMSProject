#!/bin/bash
# Program:
# 	read root file list. 
# History:
# 	2013/10/17 hchou First release
# 	2014/06/07 hchou modify
# 	2017/09/04 hchou modify
shopt -s -o nounset

if [ $# -ne 2 ]; then
    echo "illegal number of parameters."
    exit
fi

flist=$1
dataDir=$2

dpns-ls /dpm/phy.ncu.edu.tw/home/$dataDir | grep root > ./$flist
sed -i 's'',^'",root://grid71.phy.ncu.edu.tw/${dataDir}/"',g' ./$flist

echo "---- READ FILE LIST ----"
echo "FLST : " $flist
echo "PATH : " $dataDir
