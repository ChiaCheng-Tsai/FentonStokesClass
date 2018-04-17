#
# Plot results using Gnuplot
#

reset

Output=0

LW = 1
set border lw 1.4

set style line 1  linecolor rgb "black"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 2  linecolor rgb "gray"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 3  linecolor rgb "gray"  linewidth 2.000 pointtype 2 pointsize 0.8 pointinterval 5

Dir_FentonStokes = "D:/Dropbox/MyResearch/SymbolicStokesLagrange/FentonStokes/"
Dir_Fourier = Dir_FentonStokes
Dir_Stokes5 = Dir_FentonStokes
Dir_Stokes7 = Dir_FentonStokes

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
if (Output >= 0)
plot Dir_Fourier."Surface(Fourier).res" using 1:2 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Surface(Stokes5).res" using 1:2 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Surface(Stokes).res" using 1:2 title "Stokes theory 7" with linespoints ls 3;\


###############################
# Generic Velocity profiles (u)
###############################

#Xscale = 1; Yscale = 0.7

set nokey
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "u";
if (Output>0) set xlabel "$(u)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
if (Output >= 0)
plot Dir_Fourier."Flowfield(Fourier).res" using 2:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield(Stokes5).res" using 2:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield(Stokes).res" using 2:1 title "Stokes theory 7" with linespoints ls 3;\

###############################
# Generic Velocity profiles (v)
###############################

#Xscale = 1; Yscale = 0.7

set key top center reverse Left horizontal width +4
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "v";
if (Output>0) set xlabel "$(v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
if (Output >= 0)
plot Dir_Fourier."Flowfield(Fourier).res" using 3:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield(Stokes5).res" using 3:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield(Stokes).res" using 3:1 title "Stokes theory 7" with linespoints ls 3;\

unset multiplot

pause -1 "multiplot"

unset output
