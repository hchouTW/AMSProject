#!/bin/bash
# Program:
# 	read root file list. 
# History:
# 	2013/10/17 hchou First release
# 	2014/06/07 hchou modify
shopt -s -o nounset

keyWord='root'
dataDir=$1
fileList=$2

#ls $dataDir | grep $keyWord > ./$fileList
#sed -i 's'',^'",${dataDir}/"',g' ./$fileList

/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select ls $dataDir | grep $keyWord > ./$fileList
sed -i 's'',^'",root://eosams.cern.ch//${dataDir}/"',g' ./$fileList

echo "---- READ FILE LIST ----"
echo "KEY_WORD : " $keyWord
echo "PATH     : " $dataDir
echo "FILE     : " $fileList
