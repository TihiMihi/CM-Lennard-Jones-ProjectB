# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 15:08:29 2021
@author: Tihana Stefanic
s1909403
Exercise 2:
    
    An instance describes a particle in Euclidean 3D space: 
    velocity and position are [3] arrays

    Includes time integrator methods +...
"""

import math
import numpy as np
import random

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    pos: position of the particle
    vel: velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_p3d - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    """
    

    def __init__(self, label, mass, position, velocity):
        """
        Initialises a particle in 3D space

        :param label: String w/ the name of the particle
        :param mass: float, mass of the particle
        :param position: [3] float array w/ position
        :param velocity: [3] float array w/ velocity
        """
        self.label = str(label)
        self.mass = float(mass)
        self.pos = np.array(position, float)
        self.vel = np.array(velocity, float)

    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        
        :return xyz_string: string of label, x, y, z
        """

        xyz_string = str(f"{self.label} {self.pos[0]} {self.pos[1]} {self.pos[2]}")
        return xyz_string

    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        """
        KE = float(0.5 * self.mass * np.linalg.norm(self.vel)**2)
        return KE

    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        
        :return p: numpy array m*v
        """
        p_mom = np.multiply(self.mass, self.vel)
        return p_mom

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        self.vel += dt*force/self.mass


    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.pos += dt*self.vel


    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        self.pos += dt*self.vel + 0.5*dt**2*force/self.mass

    @staticmethod
    def new_p3d(file_handle):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one line per planet in the following format:
        label   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param inputFile: Readable file handle in the above format

        :return Particle3D instance
        """

        line = file_handle.readline().split()
        label = str(line[0])
        mass = float(line[1])
        position = np.array([(line[2]), (line[3]), (line[4])],float)
        velocity = np.array([float(line[5]), float(line[6]), float(line[7])])
        
        # return Particle3D instance with defined attributes label, mass, position, velocity
        
        return Particle3D(label, mass, position, velocity)

    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Computes the total kinetic energy of the system
        
        :param p3d_list: list in which each item is a P3D instance
        :return sys_ke: float of the total sum of kinetic energy of each instance
        """
        N = (len(p3d_list))
        
        sys_ke = 0
        for i in range(N):
            sys_ke += 0.5 * p3d_list[i].mass * np.linalg.norm(p3d_list[i].vel)**2
            
        return sys_ke

    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        N = len(p3d_list)
        
        total_mass = 0
        total_mom = 0
        
        for i in range(N):
            total_mass += p3d_list[i].mass
            total_mom += p3d_list[i].momentum()
        com_vel = np.divide(total_mom, total_mass)
        
        return total_mass, com_vel