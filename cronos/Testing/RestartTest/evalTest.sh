dump1=$(mktemp)
dump2=$(mktemp)

h5dump RestartTest_first_float/RestartTest_first_flt_step1.h5 > $dump1 
h5dump RestartTest_second_float/RestartTest_second_flt_step2.h5 > $dump2

nDiff=$(diff $dump1 $dump2 | wc -l)

if [ "$nDiff" -le "16" ]; then
  echo "OK"
  exit 0
else
  echo "FAILED: $nDiff"
  diff $dump1 $dump2
fi
