set terminal png background '#000000' font "Liberation 12"
set border lc rgb '#fbffff'
set border linewidth 1.5
set linetype 1 lc rgb '#86b6dc' lw 1.5
set linetype 2 lc rgb '#b3e8dc' lw 1.5
set linetype 3 lc rgb '#fff1c9' lw 1.5
set linetype 4 lc rgb '#a4f0ff' lw 1.5
set linetype 5 lc rgb '#ccf3ff' lw 1.5
set linetype 6 lc rgb '#d8fdff' lw 1.5
set linetype cycle 6
set key textcolor '#fbffff'
set title textcolor '#fbffff'
set xlabel textcolor '#fbffff'
set ylabel textcolor '#fbffff'
set output ARG1


set title "Input"
set xlabel "time"
set ylabel "amplitude"

set xrange [172000:174000]
plot "recording_1.dat" using 1 w l,\
	"recording_1.dat" using 2 w l,\
	"recording_1.dat" using 3 w l,\
	"recording_1.dat" using 4 w l,\
	"recording_1.dat" using 5 w l,\