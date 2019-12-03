from ase.io import read
from ase.calculators.dftb import Dftb


GEO_PATH = './H2O_cluster.gen'

def main():
    '''Main driver routine.'''
    system = read(GEO_PATH, format='gen')

    calc = Dftb(label='H2O_cluster', atoms=system,
                Hamiltonian_SCC='Yes',
                Hamiltonian_SCCTolerance='1.00E-010',
                Driver_='ConjugateGradient',
                Driver_MaxForceComponent='1.00E-008',
                Driver_MaxSteps=1000,
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_O='"p"',
                Hamiltonian_MaxAngularMomentum_H='"s"')

    system.set_calculator(calc)
    calc.calculate(system)

    final = read('geo_end.gen')

    forces = system.get_forces()
    energy = system.get_potential_energy()

main()

