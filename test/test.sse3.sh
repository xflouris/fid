for i in `seq 2 70`
do
  for j in `seq 2 $i`
  do
    echo $i $j
    ./tester-fid1-diag-sse3 $i $j
  done
done
