import numpy as np
import dftbplus


LIB_PATH = '/home/user/libdftbplus'


def main():
    '''Main driver routine.'''

    # coordinates of H2O, in atomic units
    coords = np.array([
        [0.000000000000000E+00, -0.188972598857892E+01,  0.000000000000000E+00],
        [0.000000000000000E+00,  0.000000000000000E+00,  0.147977639152057E+01],
        [0.000000000000000E+00,  0.000000000000000E+00, -0.147977639152057E+01]])

    # the values of extpot and extpotgrad used here were
    # taken from file: test/api/mm/testers/test_extpot.f90
    extpot = np.array([-0.025850198503435,
                       -0.005996294763958,
                       -0.022919371690684])

    extpotgrad = np.array([
        [0.035702717378527,  0.011677956375860, 0.009766745155626],
        [0.023243271928971, -0.000046945156575, 0.004850533043745],
        [0.016384005706180,  0.004608295375551, 0.005401080774962]])

    cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                              hsdpath='dftb_in.hsd',
                              logfile='log.log')

    # set geometry
    cdftb.set_geometry(coords, latvecs=None)

    # set external potential and its gradients
    cdftb.set_external_potential(extpot, extpotgrad=extpotgrad)

    # get number of atoms
    natoms = cdftb.get_nr_atoms()

    # calculate energy, forces and Gross charges
    merminen = cdftb.get_energy()
    gradients = cdftb.get_gradients()
    grosschgs = cdftb.get_gross_charges()

    # finalize DFTB+ and clean up
    cdftb.close()

    # print obtained nr. of atoms
    print('(H2O) Obtained nr. of atoms: {:d}'.format(natoms))

    # print obtained mermin free energy
    print('(H2O) Obtained Mermin-energy: {:15.10f}'.format(merminen))

    # print obtained gradients
    print('(H2O) Obtained gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[0]))
    print('(H2O) Obtained gradient of atom 2: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[1]))
    print('(H2O) Obtained gradient of atom 3: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[2]))

    # print obtained Gross charges
    print('(H2O) Obtained Gross charges: {:15.10f} {:15.10f}'
          .format(*grosschgs))


if __name__ == "__main__":
    main()
