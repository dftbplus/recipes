from ase.io import read, write
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb


NIMAGES = 13


def main():
    '''Main driver routine.'''

    initial = read('NH3_initial.traj')
    final = read('NH3_final.traj')

    images = [initial]
    images += [initial.copy() for ii in range(NIMAGES)]
    images += [final]

    neb = NEB(images)
    neb.interpolate()

    opt = BFGS(neb, trajectory='i2f.traj')

    calcs = [Dftb(label='NH3_inversion',
                  Hamiltonian_SCC='Yes',
                  Hamiltonian_SCCTolerance='1.00E-06',
                  Hamiltonian_MaxAngularMomentum_N='"p"',
                  Hamiltonian_MaxAngularMomentum_H='"s"')
             for ii in range(NIMAGES)]

    for ii, calc in enumerate(calcs):
        images[ii + 1].set_calculator(calc)

    opt.run(fmax=1.00E-02)


if __name__ == "__main__":
    main()
