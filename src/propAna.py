from MDAnalysis import Universe 
import MDAnalysis as mda
import sys 
import numpy as np
from tqdm import tqdm 
import argparse

class EnergyFlowAnalyzer():

    def __init__(self):
        pass
    
    def load_traj(self, ref, trr):
        self.traj = Universe(ref, trr)
        return self
    
    def kinetic_energy(self, mass, v, fdeg=1.5):
        """
        fdeg: 1.5 = 3/2
        """
        return 3.0/2.0 * mass * np.dot(v,v)

    def get_velocities_of_each_residue(self, selection="protein"):
        sele = self.traj.select_atoms(selection)
        if len(sele.atoms) != len(sele.velocities):
            sys.exit("Error: The number of atoms is not equal to the number of velocies.")

        with mda.Writer(f"traj.pdb", 'w') as w:   
            for ts in tqdm(self.traj.trajectory[1:]):

                for i, velocity in enumerate(sele.velocities):
                    atom = sele.atoms[i]
                    K = self.kinetic_energy(atom.mass, velocity)

                    # three-value kinetic energy is assigned with each atom
                    if K <=5000.0:
                        atom.tempfactor = 0.0

                    elif K > 5000.0 and K <= 7000.0:
                        atom.tempfactor = 0.5

                    else:
                        atom.tempfactor = 1.0

                w.write(sele.atoms)

        return self
    
def main():
    p = argparse.ArgumentParser()
    p.add_argument("-f", "--trr", required=True)
    p.add_argument("-s", "--ref", required=True)
    args = p.parse_args()
    trr = args.trr
    ref = args.ref

    EFA  = EnergyFlowAnalyzer()
    EFA.load_traj(ref, trr)
    EFA.get_velocities_of_each_residue()

if __name__ == "__main__":
    main()
