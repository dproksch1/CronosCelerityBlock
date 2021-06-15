#!/bin/bash

# new version of batchXdmfGen to create single XMF file
# (w) JK 11jul12 -- last change 07jul19

if [ $# -eq 2 ]; then
    poub=$1
    pname=$2
    cur_dir=$PWD
    cd $poub/$pname\_float
    echo $PWD
    pre=$pname\_flt_step
    suf=.h5
    itimes=$(ls $pre*$suf | grep -v _coord | sed -e 's/'$pre'\(.*\)'$suf'/\1/')
    cd $cur_dir
    echo $itimes
    ~/src/CronosCode/cronos/DataReader/Linux-amd64/XdmfGen $poub $pname $itimes
else
    echo " Usage: batch_XdmfGen <poub> <pname>"
fi
