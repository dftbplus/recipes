from ase.io import read
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb


GEO_PATH = './H2O_cluster.gen'

def main():
    '''Main driver routine.'''
    system = read(GEO_PATH, format='gen')

    system.set_calculator(Dftb(label='H2O_cluster', atoms=system,
                               Hamiltonian_SCC='Yes',
                               Hamiltonian_SCCTolerance=1.00E-010,
                               Hamiltonian_MaxAngularMomentum_='',
                               Hamiltonian_MaxAngularMomentum_O='"p"',
                               Hamiltonian_MaxAngularMomentum_H='"s"'))

    opt = BFGS(system, trajectory='opt.traj', logfile='opt.log')
    opt.run(fmax=1.00E-008)

    forces = system.get_forces()
    energy = system.get_potential_energy()

main()

