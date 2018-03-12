#!/bin/bash
output=$1
if [[ ${output} == "" ]]; then exit; fi

dpmhost=${DPM_HOME}/${output}
filelist=$(dpns-ls -lh ${dpmhost} | awk '{print $5 " " $9}' | grep "0 " | awk '{print $2}')
for file in $filelist
do
    path=${dpmhost}/${file}
    echo "dpns-rm $path"
    dpns-rm $path
done
