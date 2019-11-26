#!/bin/bash
# Local Job Command
# ljsearch
# ljkill
# ljcheck

function ljsearch {
    if [ $# == 1 ]; then
        ps -U $USER -o ppid,user,s,pcpu,pmem,cputime,cmd | head -n 1
        ps -U $USER -o ppid,user,s,pcpu,pmem,cputime,cmd | grep ${1} | grep -v grep
    else
        ps -U $USER -o ppid,user,s,pcpu,pmem,cputime,cmd
    fi
    echo ""
}

function ljkill {
    if [ $# == 1 ]; then
        ps -U $USER -o pid,cmd | grep ${1} | grep -v grep | awk '{print $1}' | xargs kill
    else
        echo -e "Error: No Keyword."
    fi
    echo ""
}

function ljcheck {
    if [ $# == 1 ]; then
        date_beg=`date`
        while true
        do
            clear
            jobs_num=$((`ljsearch ${1} | wc -l`-2))
            echo -e "DATE BEGIN ${date_beg} NOW `date`"
            echo -e "NJOBS ${jobs_num}\n"
            ljsearch ${1}
            if (( ${jobs_num} == 0 )); then
                echo -e "SUCCESS."
                break
            else
                sleep 30
            fi
        done
    else
        echo -e "Error: No Keyword."
    fi
    echo ""
}
