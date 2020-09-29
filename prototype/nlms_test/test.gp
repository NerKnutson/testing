set terminal png
set output 'test_output.png'
plot "test_output.dat" using 1,\
  "" u 2,\
  "" u 3,\
  "" u 4,\
  "" u 5
