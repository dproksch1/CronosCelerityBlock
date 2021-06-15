#!/bin/bash

# batch-run python script to extract phase velocity and decay parameter
# and output table for further processing in gnuplot
#                                     (w) JK 25sep17

export targetfile=multiwaves-results.dat  # result text file

if [ -z "$poub" ]; then
    echo "ERROR: 'poub' var is unset or empty."
else
    projdir=$(pwd)
    cd $poub
    ls test-multiwaves-*.cat > _cat_list.txt
    cd $projdir
    mv $poub/_cat_list.txt .
    sed -e 's/.cat//' _cat_list.txt > _proj_list.txt
    cat _proj_list.txt
    echo -n "Process results from these simulations? Hit # to confirm. "
    read conf
    if [ "$conf" = "#" ]; then
	echo "## k      v_phase    damping " >  $targetfile
        echo "##---------------------------" >> $targetfile
	cat _proj_list.txt | while read ksim  # loop over entries
	do
	    python analysis/wavefit.py $ksim | grep "##" > _extr.dat
	    cat _extr.dat
	    sed -e 's/##//' -e 's/k=//' -e 's/v_phase=//' \
		-e 's/damping=//' _extr.dat >> $targetfile
	done
	echo "--> data written to file \"$targetfile\"."
    else
	echo "nothing done."
    fi
    rm _cat_list.txt _proj_list.txt _extr.dat
fi
