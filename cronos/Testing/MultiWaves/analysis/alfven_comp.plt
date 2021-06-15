# Generate plots of Fig.12 in Cronos ref. (Kissmann+2018)
# requires: files "damp_040-multi.dat", "damp_040-naive.dat", and
            "damp_alfven.dat" in the same folder
# output: Postscript file "test-compare-damping.ps"
#
# Usage:  gnuplot alfven_comp.plt 
# Last modified: 13feb2018 <jk>

#set encoding iso_8859_1
set term postscript enhanced solid color
set size 0.7
set termoption dash

set style line 1 lt 3 lc rgb "blue" lw 2
set style line 2 lt 2 lc rgb "red" lw 2
set style line 3 lt 1 lc rgb "black" lw 2

set output "test-compare-damping.ps"

set title "Damping rate (k = 4.0)"
set key left bottom
#set key reverse
set xlabel "time"
set ylabel "B_z amplitude"

plot [0:2][0:0.06] \
    "damp_040-multi.dat" u 1:3 ls 1 t "Two-fluid proper init", \
    "damp_040-naive.dat" u 1:3 ls 2 t "Two-fluid 'naive' init" , \
    "damp_alfven.dat"    u 1:3 ls 3 t "Single-fluid Alfven",

set terminal x11
