#!/bin/bash
if [ $# -ne 1 ]; then
    echo "illegal number of parameters."
    exit
fi

output=$1
dpmhost=${DPM_HOME}/${output}
filelist=$(dpns-ls -lh ${dpmhost} | awk '{print $5 " " $9}' | grep "0 " | awk '{print $2}')
for file in $filelist
do
    path=${dpmhost}/${file}
    echo "dpns-rm $path"
    dpns-rm $path
done
