#!/bin/env python                                                                                                                                                             
# Author: Lin (lin1209@purdue.edu)
import sys
sys.path.append('../../')
from frame_generator import frame_generator

def main(argv):
   
      
   # print out all atoms' x coordinate and time step
   for atom,timestep,box in frame_generator('0.nvt.lammpstrj',end=1):
         print(atom.x)
         print(timestep)



if  __name__ == '__main__': 
    main(sys.argv[1:])
