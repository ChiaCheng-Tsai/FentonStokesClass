#
# Plot results using Gnuplot
#

reset

Output=0

LW = 1
set border lw 1.4

set style line 1  linecolor rgb "black"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 2  linecolor rgb "gray"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 3  linecolor rgb "green"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0

Dir_Fourier = "D:/research/codes/Fourier/"
Dir_Stokes5 = Dir_Fourier."Stokes/"
Dir_Stokes7 = Dir_Fourier."Stokes7/"

set multiplot layout 3,1

#########################
# Generic Surface profile
#########################

set nokey
set format x "%3.0f"
set format y "%4.1f"
if (Output==0) set xlabel "x/h"; set ylabel "y/h"
if (Output>0) set xlabel "$x/d$" offset 0,-0.5; set ylabel "$\\dfrac{y}{d}$" offset -1
set xtics nomirror
set ytics 0, 0.5, 2
#show ytics
set yrange [0:1.6]
set autoscale x
File = Dir_Stokes7."Surface"
if (Output >= 0)
plot Dir_Fourier."Surface.res" using 1:2 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Surface.res" using 1:2 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Surface.res" using 1:2 title "Stokes theory 7" with lines ls 3;\


###############################
# Generic Velocity profiles (u)
###############################

#Xscale = 1; Yscale = 0.7

set nokey
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Velocity u (dimensionless w.r.t. g & h)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
File = Dir_Stokes7."u"
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 2:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 2:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 2:1 title "Stokes theory 7" with lines ls 3;\

###############################
# Generic Velocity profiles (v)
###############################

#Xscale = 1; Yscale = 0.7

set key top center reverse Left horizontal width +4
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Velocity v (dimensionless w.r.t. g & h)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
File = Dir_Stokes7."v"
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 3:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 3:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 3:1 title "Stokes theory 7" with lines ls 3;\

unset multiplot

pause -1 "multiplot"

unset output
