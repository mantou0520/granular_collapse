# granular_collapse
code for testing the collapse of granular columns

In this practice, we use a library called Mechsys (http://mechsys.nongnu.org/)

column.cpp contains everything about the setup of the granular packing, the shape of the particles and shape 
of the granular columns.

column.inp is the input file

Once built, use the following command to run
./column column 48
where 48 is the number of core we are going to use.

Also, we could change the shape of the cross-section of the granular column to either circular, or right triangle, or isoscele triangle. This can help us investigate the influence of the initial condition on the dynamics of the granular collapsing.
