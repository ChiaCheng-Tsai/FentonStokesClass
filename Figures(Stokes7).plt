#
# Plot results using Gnuplot
#

reset

Output = 0

LW = 1
set border lw 1.4

Xscale = 1; Yscale = 0.4
#set style line 1 linetype 1 linecolor rgb "blue" lw LW
#set style line 2 linetype 1 linecolor rgb "red" lw LW
#set style line 3 linetype 2 linecolor rgb "green" lw LW
#set style line 4 linetype 2 linecolor rgb "magenta" lw LW
#set style line 6 linetype 6 linecolor rgb "black" ps 0.6 lw 2.5

set style line 1  linecolor rgb "black"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 2  linecolor rgb "gray"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0
set style line 3  linecolor rgb "gray"  linewidth 2.000 pointtype 2 pointsize 2 pointinterval 0

Dir_Fourier = "D:/Dropbox/MyResearch/Fourier/"
Dir_Stokes5 = Dir_Fourier."Stokes/"
Dir_Stokes7 = Dir_Fourier."Stokes7/"

#########################
# Generic Surface profile
#########################

set key bottom center reverse Left horizontal width +8
set format x "%3.0f"
set format y "%4.1f"
if (Output==0) set xlabel "x/d"; set ylabel "y/d"
if (Output>0) set xlabel "$x/d$" offset 0,-0.5; set ylabel "$\\dfrac{y}{d}$" offset -1
set xtics nomirror
set ytics 0, 0.5, 2
#show ytics
set yrange [0:1.6]
set autoscale x
File = Dir_Stokes7."Surface"

load Dir_Fourier.'SetOutput.plt'
if (Output >= 0)
plot Dir_Fourier."Surface.res" using 1:2 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Surface.res" using 1:2 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Surface.res" using 1:2 title "Stokes theory 7" with lines ls 3;\
     pause -1 "Surface"

###############################
# Generic Velocity profiles (u)
###############################

Xscale = 1; Yscale = 0.7

set key top center reverse Left horizontal width +4
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Velocity u (dimensionless w.r.t. g & d)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
File = Dir_Stokes7."u"
load Dir_Fourier.'SetOutput.plt'
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 2:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 2:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 2:1 title "Stokes theory 7" with lines ls 3;\
pause -1 "u"

###############################
# Generic Velocity profiles (v)
###############################

Xscale = 1; Yscale = 0.7

set key top center reverse Left horizontal width +4
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Velocity v (dimensionless w.r.t. g & d)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
File = Dir_Stokes7."v"
load Dir_Fourier.'SetOutput.plt'
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 3:1 title "Fourier approximation" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 3:1 title "Stokes theory 5" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 3:1 title "Stokes theory 7" with lines ls 3;\
pause -1 "v"

############################
# Generic Velocity profiles
############################

Xscale = 1; Yscale = 0.7

set key top center reverse Left horizontal width +4
set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Velocity (dimensionless w.r.t. g & d)";
if (Output>0) set xlabel "$(u,v)/\\sqrt{g d}$" offset 0,-0.5;
set title ""
File = Dir_Stokes7."uv"
load Dir_Fourier.'SetOutput.plt'
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 2:1 title "Fourier - u" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 2:1 title "Stokes 5 - u" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 2:1 title "Stokes 7 - u" with lines ls 3,\
     Dir_Fourier."Flowfield.res" using 3:1 title "Fourier - v" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 3:1 title "Stokes 5 - v" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 3:1 title "Stokes 7 - v" with lines ls 3;\
pause -1 "uv"

############################
# d(Velocity)/dt profiles
############################

set format x "%5.2f"
set format y "%4.1f"
if (Output==0) set xlabel "Partial du/dt and dv/dt (dimensionless w.r.t. g & d)";
if (Output>0) set xlabel "${\\partial u}/{\\partial t}$ and ${\\partial v}/{\\partial t}$  (dimensionless w.r.t. $g$ \\& $d$)" offset 0,-0.5;
set title "Time derivatives of velocity"
File = Dir_Stokes7."utvt"
load Dir_Fourier.'SetOutput.plt'
if (Output >= 0)
plot Dir_Fourier."Flowfield.res" using 5:1 title "Fourier - H" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 5:1 title "Stokes 5 - H" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 5:1 title "Stokes 7 - H" with lines ls 3,\
     Dir_Fourier."Flowfield.res" using 6:1 title "Fourier - V" with lines ls 1,\
     Dir_Stokes5."Flowfield.res" using 6:1 title "Stokes 5 - V" with lines ls 2,\
     Dir_Stokes7."Flowfield.res" using 6:1 title "Stokes 7 - V" with lines ls 3;\
     pause -1 "utvt"

unset output
