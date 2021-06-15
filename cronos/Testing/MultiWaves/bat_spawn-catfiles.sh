#!/bin/bash

# generate sequence of cat files, one for each value of k
# xe changed such that box size is L/pi = 2/k (one cycle)
# and time only covers one period.
# note 'scale_x = pi' entry!
#                                     (w) JK 22sep17

if [ -z "$poub" ]; then
    echo "ERROR: 'poub' var is unset or empty."
else            ## seq: start-stride-end, padded with zeroes
    seq -w 5 5 100 > _kvals10.txt   ## -> 005, 010, ..,, 100
    cat _kvals10.txt | while read num  # loop over entries
    do
	file=test-multiwaves-$num.cat  # target cat file name
	cp test-multiwaves-STENCIL.cat _temp1.cat
	frac=$(echo "scale=1;  $num/10." | bc) # divide by 10
	xmax=$(echo "scale=8;  2./$frac" | bc) # compute xe = 2/k
	tend=$(echo "scale=8; 10./$frac" | bc) # compute t_end = 10/k
	dtfl=$(echo "scale=8; 0.2/$frac" | bc) # compute dt_float (~50 frames)
	echo 'writing "'$poub/$file'" (k='$frac', xe='$xmax', t_end='$tend')'
	sed -e 's/<X_END>/'$xmax'/' _temp1.cat > _temp2.cat   # insert xe
	sed -e 's/<T_END>/'$tend'/' _temp2.cat > _temp1.cat   # insert t_end
	sed -e 's/<DT_FL>/'$dtfl'/' _temp1.cat > $poub/$file  # insert dt_float
    done
    rm _kvals10.txt _temp1.cat _temp2.cat
    echo "...done."
fi
