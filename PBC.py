# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 15:11:35 2022
@author: Tihana
s1909403
Exercise 1: PBC
"""
import numpy as np
import random

def pbc(x,l):
    """ PBC: Return the image/position of x_0 in the original box where 0 < x < l:
        simulation in the original box is determined by the remainder of (vector/boxlength) from the origin
    """
    x_0 = np.mod(x,l)
    return x_0

def mic(x,l):
    """ MIC: Return the image/position of the nearest repeating x_0 particle to the origin:
        Nearest neighbour has coordinates closest to (0,0,0) which is equivalent to original x_0 from its closest cube coordinate n(l,l,l)
        the minimum image convention is determined by comparing the remainder of vector x / boxlength l to the boxlength
    """
    x_0 = np.mod((x+(l/2)),l) - l/2
    return x_0

def random_v():
    """ Return randomized values for vector x and boxlength l """
    l = random.randint(1, 10)
    x = np.random.uniform(low=0, high=20, size=3)
    """ Define variables x and l """
    return x,l
x,l = random_v()

def main(x,l):
    """ Test the functions for randomly generated vector x and boxlength l """
    print(f"\nFor randomly generated input of variables: \nBox length = {l} \nVector x = {x}")
    print(f'x_0 in periodic boundary conditions, image in the original box 0 < x_0 < l: \nPBC = {pbc(x, l)}')
    print(f"x_0 in the minimum image convention, the nearest neighbour to the origin: \nMIC = {mic(x,l)}")
    
    """ Test the functions for given values of vector x and boxlength l """
    x = np.array([5, 7, 7])
    l = 6
    print(f"\nFor vector x = {x} and box length l = {l}:")
    print(f"x_0 in periodic boundary conditions, image in the original box 0 < x_0 < l: \nPBC = {pbc(x, l)}")
    print(f"x_0 in the minimum image convention, the nearest neighbour to the origin: \nMIC = {mic(x,l)}\n")

    x = np.array([12, 13, 14])
    l = 8
    print(f"For vector x = {x} and box length l = {l}:")
    print(f"x_0 in periodic boundary conditions, image in the original box 0 < x_0 < l: \nPBC = {pbc(x, l)}")
    print(f"x_0 in the minimum image convention, the nearest neighbour to the origin: \nMIC = {mic(x,l)}")

    return x, l, pbc(x,l), mic(x,l)

# Execute main method, but only if it is invoked directly
if __name__ == "__main__":
    main(x,l)