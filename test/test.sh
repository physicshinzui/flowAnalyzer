#!/bin/bash 

python ../src/propAna.py -f ../mod.trr -s ../ref.pdb
pymol ../src/load.py -- traj.pdb
