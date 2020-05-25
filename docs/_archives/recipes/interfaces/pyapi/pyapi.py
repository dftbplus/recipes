import numpy as np
import dftbplus


LIB_PATH = '/home/user/libdftbplus'

# DFTB+ conversion factors
# (according to prog/dftb+/lib_common/constants.F90)
BOHR__AA = 0.529177249
AA__BOHR = 1 / BOHR__AA


def main():
    '''Main driver routine.'''

    # coordinates of TiO2, in atomic units
    coords = np.array([
        [-0.016726922839251,  0.016725329441158, -0.000003204152532],
        [-0.016726505918979,  1.920201169305565, -7.297102897292027],
        [ 0.017412997824265, -0.024318617967798,  2.005339137853385],
        [ 1.920770753428742, -0.024319922392223, -4.437737763954652],
        [ 0.024319174400169, -0.017404302527510, -2.005347277168561],
        [ 0.024317270342179,  1.886164739806594, -5.291732430733527]])
    coords *= AA__BOHR

    # lattice vectors of TiO2, in atomic units
    latvecs = np.array([
        [-1.903471721000000,  1.903471721000000,  4.864738245000000],
        [ 1.903471721000000, -1.903471721000000,  4.864738245000000],
        [ 1.903471721000000,  1.903471721000000, -4.864738245000000]])
    latvecs *= AA__BOHR

    cdftb = dftbplus.DftbPlus(libpath=LIB_PATH,
                              hsdpath='dftb_in.hsd',
                              logfile='TiO2.log')

    # set geometry
    cdftb.set_geometry(coords, latvecs=latvecs)

    # get number of atoms
    natoms = cdftb.get_nr_atoms()

    # calculate energy, gradients and Gross charges
    merminen = cdftb.get_energy()
    gradients = cdftb.get_gradients()
    grosschgs = cdftb.get_gross_charges()

    # finalize DFTB+ and clean up
    cdftb.close()

    # print obtained nr. of atoms
    print('(TiO2) Obtained nr. of atoms: {:d}'.format(natoms))

    # print obtained mermin free energy
    print('(TiO2) Obtained Mermin-energy: {:15.10f}'.format(merminen))

    # print obtained gradients
    print('(TiO2) Obtained gradient of atom 1: {:15.10f} {:15.10f} {:15.10f}'
          .format(*gradients[0]))

    # print obtained Gross charges
    print('''(TiO2) Obtained Gross charges: {:15.10f} {:15.10f}
          {:15.10f} {:15.10f} {:15.10f} {:15.10f}'''.format(*grosschgs))


if __name__ == "__main__":
    main()
