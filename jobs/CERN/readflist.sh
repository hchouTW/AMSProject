#!/bin/bash
# Program:
# 	read root file list. 
# History:
# 	2013/10/17 hchou First release
# 	2014/06/07 hchou modify
# 	2017/09/04 hchou modify
shopt -s -o nounset

if [ $# -lt 2 ] && [ $# -ne 3 ]; then
    echo "illegal number of parameters."
    exit
fi

head=root://eosams.cern.ch/
if [ $# -eq 3 ]; then
    head=$3
fi

flist=$1
dataDir=$2

ls $dataDir | grep root > ./$flist
sed -i 's'',^'",${head}${dataDir}/"',g' ./$flist

echo "---- READ FILE LIST ----"
echo "FLIST : " $flist
echo "PATH  : " $dataDir
