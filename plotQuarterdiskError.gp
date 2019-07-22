set terminal postscript enhanced color
set output "quarterdiskError.eps"
set title "quarterdisk function: sampling error" font ",30"
set xlabel "samples" font ",30"
set ylabel "error" font ",30"
set tics font ",25"
set key spacing 1.5
set logscale x
set logscale y

plot [16:1024] [0.001:0.1]\
 "errors_quarterdisk_random.data" using 1:2 smooth unique lw 5.0 title "random",\
 "errors_quarterdisk_bestcand.data" using 1:2 smooth unique lw 5.0 title "best cand",\
 "errors_quarterdisk_irrational_rot.data" using 1:2 smooth unique lw 5.0 title "irrational rot",\
 "errors_quarterdisk_halton_base23_owen.data" using 1:2 smooth unique lw 5.0 title "Halton Owen",\
 "errors_quarterdisk_sobol_dim01_owen.data" using 1:2 smooth unique lw 5.0 title "Sobol Owen",\
 "errors_quarterdisk_pmj02.data" using 1:2 smooth unique lw 5.0 lt rgb "black" title "pmj02",\
 0.35*x**(-0.5) lw 3.0 lt rgb "gray" title "N^{-0.5}",\
 0.4*x**(-0.75) lw 3.0 lt rgb "gray" title "N^{-0.75}"
