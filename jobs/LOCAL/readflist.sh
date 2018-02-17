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

keyWord='root'
flist=$1
dataDir=$2

ls $dataDir | grep $keyWord > ./$flist
sed -i 's'',^'",${dataDir}/"',g' ./$flist

echo "---- READ FILE LIST ----"
echo "KEY_WORD : " $keyWord
echo "FILE     : " $flist
echo "PATH     : " $dataDir
