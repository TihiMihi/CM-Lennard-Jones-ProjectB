# CM-Lennard-Jones-ProjectB
Lennard Jones simulation

Created on Tue Jan 25 14:15:28 2022

@author: Tihana Stefanic

## The Lennard Jones Potential on Argon particles 

Project B: Lennard Jones simulation of Argon
particles moving in a double well potential using
a velocity Verlet time integration.

Produces trajectory file of N particles interacting
via the Lennard Jones potential at a given density
and temperature.

Produces plots of the energies of the system
and the mean distribution as functions of time
and the average radial distribution function. 

### Instructions

*bold*
**italic**
***bold and italic***
'this is code'

Main asks for physical state user input of
particle system as shortcut to define **number
of particles**, **density**, and **temperature**

In the case of no user input or input which
does not correspond to 'solid' 'liquid' or 'gas',
main asks for manual input of system parameters

Main asks for **timestep** and **number of steps**
user input

Code automated to work with no user input for
the physical state or system parameters with
default values for a solid

Main asks for trajectory file name user input
which is also automated to name file as
'output.xyz' in the case of no user input