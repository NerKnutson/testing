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


set title "test_shot"
set xlabel "time"
set ylabel "amplitude"

plot "test_shot.dat" using 1 w l,\
	"test_shot.dat" using 2 w l,\
	"test_shot.dat" using 3 w l,\
	"test_shot.dat" using 4 w l,\
	"test_shot.dat" using 5 w l,\
	"test_shot.dat" using 6 w l,\
	"test_shot.dat" using 7 w l,\
	"test_shot.dat" using 8 w l,\
	"test_shot.dat" using 9 w l,\
	"test_shot.dat" using 10 w l,\
	"test_shot.dat" using 11 w l,\
	"test_shot.dat" using 12 w l,\
	"test_shot.dat" using 13 w l,\
	"test_shot.dat" using 14 w l,\
	"test_shot.dat" using 15 w l,\
	"test_shot.dat" using 16 w l,\
	"test_shot.dat" using 17 w l,\
	"test_shot.dat" using 18 w l,\
	"test_shot.dat" using 19 w l,\
