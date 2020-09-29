set terminal png
set output 'output.png'
plot "output.dat" using 1,\
  "" u 2,\
  "" u 3,\
  "" u 4,\
  "" u 5
