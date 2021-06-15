#!/bin/bash

# batch-run all simulations with previously created cat files
# output is suppressed
# Created:       22sep17 jk
# Last modified: 12feb18 jk (use more than one proc)

# Usage: ./bat_run-all-sims.sh <number-of-procs-to-use>

if [ -z "$poub" ]; then
    echo "ERROR: 'poub' var is unset or empty."
else
    if [ $# -eq 0 ]; then
  	let np=0
	echo 'No arg => will use just one proc.'
    else
	let nproc=$1
	echo '-> will use '$nproc' procs.'
    fi
    projdir=$(pwd)
    cd $poub
    ls test-multiwaves-*.cat > _cat_list.txt
    cd $projdir
    mv $poub/_cat_list.txt .
    sed -e 's/.cat//' _cat_list.txt > _proj_list.txt
    cat _proj_list.txt
    echo -n "Execute these simulations, removing old data? Hit # to confirm. "
    read conf
    if [ "$conf" = "#" ]; then
	make -j2
	let sim_count=0
	cat _proj_list.txt | while read ksim  # loop over entries
        do
	    rm -rf $(echo $poub/$ksim'_*')
	    let step=$(expr 1 + $sim_count % $nproc)  # in [0, 1, ..., nproc-1]
	    let sim_count=sim_count+1
	    let batch=$(expr 1 + $(expr $(expr $sim_count - $step) / $nproc))
	    echo "running" $ksim "(batch #"$batch"," $step "of" $nproc")..."
	    if [ $step -ne $nproc ]; then  # -> not last in group
		Linux-amd64/proj $poub $ksim 1> /dev/null &
	    else                           # -> last job in group
		Linux-amd64/proj $poub $ksim 1> /dev/null
	    fi
     	    		    done
	echo "ready."
    else
	echo "nothing done."
    fi
    rm _cat_list.txt _proj_list.txt
fi
