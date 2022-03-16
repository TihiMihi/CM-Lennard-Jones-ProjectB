# CM-Lennard-Jones-ProjectB
Lennard Jones simulation

Created on Tue Jan 25 14:15:28 2022

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

Main asks for a base file name user input to make
the output file names which is automated to construct 
file names using 'output' in case of no user input

### Format

#### Input parameters:
Number of particles
Density
Temperature
Timestep
Number of steps
Output files base name

#### Output files:
***Trajectory file***: label and positions of
particle system at a time instance

***Energy file***: kinetic, potential, and
total energy of particle system at a time instance

***MSD file***: mean squared displacement of
particle system at a time instance from time zero

***RDF file***: radial distribution function
averaged across all timesteps in simulation

### Units

Distances in units of sigma, where sigma is
the size of the particle

Mass in units such that the mass of a particle is 1

Energy in units of epsilon, where epsilon
is the dispersion energy, also quantifies
the depth of the potential well

Temperature in units of epsilon i.e. in terms of
the Boltzman constant

Time in units of the ratio of unitary mass
to epsilon all square rooted times sigma