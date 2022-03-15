# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 14:15:28 2022

@author: Tihana Stefanic

<h1> The Lennard Jones Potential on Argon particles </h1>

<h2> Project B: Lennard Jones simulation of Argon
particles using velocity Verlet time integration 
of a particle moving in a double well potential. </h2>

<h3> Produces trajectory file of N particles interacting
via the Lennard Jones potential at given density
and temperature. </h3>

<h4> Produces plots of the energies of the system
and the mean distribution function as functions
of time and the average radian distribution function. </h4>

"""

import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D
from PBC import mic, pbc
from mdutilities_20210131 import set_initial_positions
from mdutilities_20210131 import set_initial_velocities

def mic_separation(particle_list, cell_size):
    """
    Minimum Image Convention separation of particles

    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    cell_size : length of
        simulation cell

    Returns
    -------
    separation : array of
        particle separations

    """
    
    N = len(particle_list)
    separation = np.zeros([N,N,3])
    
    for j in range(N):
        for i in range(j):
            separation[i,j] = mic((particle_list[i].pos - particle_list[j].pos), cell_size)
            separation[j,i] = -separation[i,j]

    return separation

def new_force(particle_list, separation):
    """
    Method to return force between particles
    interacting via Lennard Jones potential
    
    Force is given by:
    F(r) = 48 * (1/r^14 - 1/2r^8)(r2-r1)
    with r1 and r2 being particle positions
    where r = |r2 - r1|

    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    separation : array of
        particle separations
    
    Returns
    -------
    force_new : array of Lennard Jones
        force for N particles
        
    """
    
    N = len(particle_list)
    force_new = np.zeros((N,3))
    
    r = separation
    r_mod = np.linalg.norm(separation, axis=-1)
    
    for j in range(N):
        for i in range(j):
            if r_mod[i,j] <= 3.5: # cut-off radius
                fij = (48*(r_mod[i,j]**(-14) - 0.5*r_mod[i,j]**(-8))*r[i,j])
                force_new[i] += fij
                force_new[j] -= fij
            
    return force_new

def new_potential(particle_list, separation):
    """
    Method to return potential energy between
    particles interacting via Lennard Jones potential
    
    Potential is given by:
    V(r) = 4 * ((1/r^12 - 1/r^6)
    with r1 and r2 being particle positions
    where r = |r2 - r1|

    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    separation : array of
        particle separations

    Returns
    -------
    potential_new : float of Lennard
        Jones Potential energy

    """
    
    N = len(particle_list)
    potential_new = 0
    
    #r_mod = np.linalg.norm(separation, axis=-1)
    r_mod = np.linalg.norm(separation, axis = 2)
    
    for j in range(N):
        for i in range(j):
            if i != j:
                potential_new += 4 * (1/r_mod[i,j]**12 - 1/r_mod[i,j]**6)
            
    return potential_new

def xyz_trajectory(particle_list, file_handle):
    """
    Method to output label and positions
    of particle3D instances in list
    :param file_handle: file handle

    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    file_handle : file handle

    Returns
    -------
    None.
    """
    
    file_handle.write(str(len(particle_list)))
    file_handle.write(f"\nComment line\n")
    
    for particle in particle_list:    
        file_handle.write(str(particle) +"\n")

def total_energies(particle_list, cell_size, separation, file_handle):
    """
    Method to kinetic, potential,
    and total energy of system
    
    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    cell_size : length of
        simulation cell
    separation : array of
        particle separations
    file_handle : file handle

    Returns
    -------
    KE : kinetic energy
        float
    PE : potential energy
        float
    total_energy : total
        energy float

    """

    KE = 0
    PE = 0
    total_energy = KE + PE
    
    for i in range(len(particle_list)):
        KE += particle_list[i].kinetic_e()
    
    PE += new_potential(particle_list, separation)
    total_energy += (KE + PE)
    
    file_handle.write(f"\nKinetic energy: {KE}\nPotential energy: {PE}\nTotal energy:{total_energy}\n")
    
    return KE, PE, total_energy

def MSD(particle_list, cell_size, initial_pos_list, time_list, iterator, file_handle):
    """
    Parameters
    ----------
    particle_list : list of 
        Particle3D instances
    cell_size : length of
        simulation cell
    initial_pos_list : list of
        initial particle positions
    time_list : list of
        time increments
    iterator : time
        iterator
    file_handle : file handle
    
    Returns
    -------
    mean square displacement
        array
    
    """
    N = len(particle_list)
    msd = 0
    time_pos_list = np.zeros((N,3))
    
    for i in range(N):
        time_pos_list[i] = particle_list[i].pos
    
    if len(time_list) == iterator: # time == time_list[iterator]
        
        delta_pos = mic(time_pos_list - initial_pos_list, cell_size)
        msd = np.sum(delta_pos**2) / N
            
        file_handle.write(f"{time_list[-1]} {msd}\n")
    
        return msd # return msd for only selection of time intervals

def RDF(N, separation_list, rho, numstep, file_handle):
    """
    
    Parameters
    ----------
    N : Number of particles.
    separation_list : list of particle
        separations overtime
    rho : density
    numstep : number of timesteps
    file_handle : file handle.
    
    Returns
    -------
    r : histogram edge
        distances
    gr : averaged weighted
        rdf histogram
    
    """
    histogram_sum = 0
    
    for separation in separation_list:
        r_mod = np.linalg.norm(separation, axis = -1)
        counts, edges = np.histogram(r_mod, bins=100, range=(0,2))
        histogram_sum += counts
    
    # time average of histograms
    avg_count = histogram_sum/numstep
    dr = edges[1]-edges[0]
    
    # list of edges excluding final edge
    # with first value annuelled
    r = edges[1:-1] + 0.5*dr
    
    # weighting
    rho_0 = 4*np.pi*rho*(r**2)*dr
    
    # average radial distribution
    gr = avg_count[1:]/(rho_0*N)
    
    file_handle.write(f"{r} \n{gr} \n")
        
    return r, gr

# Begin main code
def main():
    """
    Main function to create a list of particle3D instances of
    argon particles interacting via Lennard Jones
    
    Initial positions and velocities imported with the cell size
    Separation bewteen particles in terms of minimum image convention
    Particles set to interact via Lennard Jones potential and forces
    Update positions, separations, forces, velocities, and energies
    of particles per timestep
    
    Output trajectory results
    Output mean square displacement as function of time
    and radial distribution function as function of distance
    
    Plot total energy of the system overtime
    Plot MSD per time
    Plot RDF per distance

    Returns
    -------
    time_list : list
        list of time increment floats
    energy_list : list
        list of total energy of system in simulation
    """
    
    # System parameters with user input
    # and default values for Argon solid
    N = int(input("Number of particles:") or "32")
    print(f"N = {N}") #no of particles
    
    rho = float(input("Density in mass/(sigma^3):") or "1")
    print(f"ρ = {rho}") # density in terms of unitary mass per volume of particle
    
    temp = float(input("Temperature in terms of kB (energy):") or "0.1")
    print(f"T = {temp}") # temperature in terms of dispersion energy E
    
    
    # initialize particle list
    particle_list = [Particle3D(f"p{ii+1}", 1.0, np.zeros(3), np.zeros(3)) 
                     for ii in range(N)]
    
    
    # import cell size, initial positions and velocities
    cell_size, _ = set_initial_positions(rho, particle_list)
    set_initial_velocities(temp, particle_list)
    
    
    # print and store initial particle positions
    print("Initial particle positions:")
    
    initial_positions = np.empty(1)
    for particle in particle_list:
        np.append(initial_positions, particle.pos)
        print(particle)
    
    
    # Simulation parameters with user input
    # and default values for Argon solid
    dt = float(input("Timestep:") or "0.01")
    print(f"dt = {dt}")
    
    numstep = int(input("Number of steps:") or "1000")
    print(f"Steps = {numstep}")
    
    time = 0
    
    # Compute separations between particles
    # using MIC and store in list
    separation = mic_separation(particle_list, cell_size)
    separation_list = [separation]
    
    # Measure force
    force = new_force(particle_list, separation)
    
    # Measure energy data
    # to output to energy file
    energy_file = open("energies.dat", "w")
    energy_file.write("Energies of the system")
    # energy in terms of the dispersion energy E i.e. depth of potential well
    kinetic_energy, potential_energy, energy = total_energies(particle_list, cell_size, separation, energy_file)

    # Initialise data lists for plotting total & potential energy later
    time_list = [time]
    energy_list = [energy]
    pot_list = [potential_energy]
    
    # Output initial particle list to XYZ trajectory file
    trajectory_file = open(input("Trajectory file:") or "output.xyz", "w")
    xyz_trajectory(particle_list, trajectory_file)
    
    # Measure MSD of system as a function of time
    # to output to MSD file
    msd_file = open("MSD.dat", "w")
    msd_file.write("Time / Mean Square displacement\n")
    msd = MSD(particle_list, cell_size, initial_positions, time_list, len(time_list), msd_file)
    msd_list = [msd]
    
    
    for step in range(numstep):
        # Start the time integration loop
        
        for i in range(N):
            # Update particle position
            # using periodic boundary conditions
            particle_list[i].leap_pos2nd(dt, force[i])
            particle_list[i].pos = pbc(particle_list[i].pos, cell_size)
                        
        # Update separations between particles
        # using MIC and store in list
        separation = mic_separation(particle_list, cell_size)
        separation_list.append(separation)
        
        # Update force
        force_new = new_force(particle_list, separation)
        
        # Update velocity
        for i in range(N):
            # Update particle velocity by averaging
            # current and new forces
            particle_list[i].leap_velocity(dt, 0.5*(force[i] + force_new[i]))
    
        # Re-define force value
        force = force_new
        
        # Increase time
        time += dt
        
        # Output energy data
        kinetic_energy, potential_energy, energy = total_energies(particle_list, cell_size, separation, energy_file)
        
        # Append information to data lists
        time_list.append(time)
        energy_list.append(energy)
        pot_list.append(potential_energy)
        
        # Output particle list to XYZ trajectory file
        xyz_trajectory(particle_list, trajectory_file)
        
        # Compute MSD of system in several
        # time intervals and store in list
        msd_tlist = np.arange(1, numstep, 10)
        for i in msd_tlist:
            if i == len(time_list):
                msd = MSD(particle_list, cell_size, initial_positions, time_list, len(time_list), msd_file)
                msd_list.append(msd)
    
    # Output energy data
    kinetic_energy, potential_energy, energy = total_energies(particle_list, cell_size, separation, energy_file)
    
    # Measure radial distribution of
    # system to output to rdf file
    rdf_file = open("RDF.dat", "w")
    rdf_file.write("Time / Radial Distribution Function\n")
    r, gr = RDF(N, separation_list, rho, numstep, rdf_file)
    
    # Close files
    energy_file.close()
    trajectory_file.close()
    msd_file.close()
    rdf_file.close()
    
    # Plot potential energy
    pyplot.title('Lennard Jones: potential energy vs time')
    pyplot.xlabel('Time (sigma√(m/E)') # sigma = size of the particle, m = unitary mass, E
    pyplot.ylabel('Potential Energy (dispersion energy E)') # E = dispersion / classical binding energy
    pyplot.plot(time_list, pot_list)
    pyplot.show()
    
    # Plot total energy
    pyplot.title('Lennard Jones: total energy vs time')
    pyplot.xlabel('Time (sigma√(m/E)') # sigma = size of the particle, m = unitary mass, E
    pyplot.ylabel('Total energy (dispersion energy E)') # E = dispersion / classical binding energy
    pyplot.plot(time_list, energy_list)
    pyplot.show()
    
    # Plot Mean Square Displacement as a function of time
    pyplot.title('Lennard Jones: Mean Square Displacement vs time')
    pyplot.xlabel('Time (sigma√(m/E)') # sigma, m = unitary mass, E = dispersion energy
    pyplot.ylabel('MSD (sigma)') # sigma = size of the particle
    pyplot.plot(msd_tlist, msd_list)
    pyplot.show()
    
    # Plot Radial Distribution Function as a function of distance
    pyplot.title('Lennard Jones: Radial Distribution Function vs distance')
    pyplot.xlabel('Distance (sigma)') # sigma = size of the particle, m = unitary mass, E = dispersion energy
    pyplot.ylabel('Radial Distribution Function') # epsilon = classical binding energy # depth of potential well
    pyplot.plot(r, gr)
    pyplot.show()
    
    return time_list, energy_list

# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()