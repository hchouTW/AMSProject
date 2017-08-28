#!/bin/bash
# Program:
# 	move file to trash 
# History:
# 	2014/06/11 hchou First release
shopt -s -o nounset

mv "${@}" ~/.trash
