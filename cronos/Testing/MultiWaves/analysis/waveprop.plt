# Generate plots of Fig.11 in Cronos ref. (Kissmann+2018)
# requires: files "tab_k-omRE-omIM.dat" and "multiwaves-results.dat"
#           in parent folder
# output: (1) Postscript file "test-multiwaves-v_phase.ps"
#         (2) Postscript file "test-multiwaves-damping.ps"
#
# Usage:  gnuplot waveprop.plt 
# Last modified: 13feb2018 <jk>

#set encoding iso_8859_1
set term postscript enhanced solid color
set size 0.7
set termoption dash

set style line 1 lt 2 lc rgb "blue"  lw 2
set style line 2 lt 1 lc rgb "red"   lw 2
set style line 3 lt 1 lc rgb "black" lw 2

set key left bottom
set xlabel "k"
       
set output "test-multiwaves-v_phase.ps"
set title "Phase speed"
set ylabel "Re({/Symbol w}) / k"
plot [a=0:10][0.5:1.5] \
     "../multiwaves-results.dat" u 1:2    ls 3     t "Simulation (v_k)" , \
     "../tab_k-omRE-omIM.dat"    u 1:($2/$1) ls 2 w l t "Theory (two-fluid)" , \
      sqrt(4-a*a/16.)/2.                  ls 1     t "Theory (one-fluid)"

set key right bottom
set output "test-multiwaves-damping.ps"
set title "Damping rate"
set ylabel "- Im({/Symbol w})
plot [a=0:10][0:1.2] \
     "../multiwaves-results.dat" u 1:3  ls 3  t "Simulation ({/Symbol G}_k)", \
     "../tab_k-omRE-omIM.dat" u 1:(-$3) ls 2 w l t "Theory (two-fluid)" , \
      a*a/8.                            ls 1     t "Theory (one-fluid)"

set terminal x11
