import os
from subprocess import Popen
from ase.io import read, write
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator


NIMAGES = 13

DFTBP_PATH = 'dftb+'
HSD_INPUT = 'dftb_in.hsd'
GEO_PATH = 'NH3_initial.gen'


def main():
    '''Main driver routine.'''

    system = read(GEO_PATH, format='gen')
    write('dummy.gen', system, format='gen')

    initial = read('NH3_initial.traj')
    final = read('NH3_final.traj')

    images = [initial]
    images += [initial.copy() for ii in range(NIMAGES)]
    images += [final]

    neb = NEB(images)
    neb.interpolate()

    opt = BFGS(neb, trajectory='i2f.traj')

    socketids = range(1, NIMAGES + 1)
    wdirs = ['_calc/image_{:d}'.format(socket) for socket in socketids]
    unixsockets = ['dftbplus_{:d}'.format(socket) for socket in socketids]

    write_modhsd(socketids)

    calcs = [SocketIOCalculator(log='socket.log', unixsocket=unixsocket)
             for unixsocket in unixsockets]

    for ii, calc in enumerate(calcs):
        images[ii + 1].set_calculator(calc)

    for cwd in wdirs:
        Popen(DFTBP_PATH, cwd=cwd)

    opt.run(fmax=1.00E-02)

    for calc in calcs:
        calc.close()


def write_modhsd(socketids):
    '''Writes input files based on 'dftb_in.hsd' with modified unixsockets.

    Args:
        socketids (1darray): contains unixsocket identifier
    '''

    for socket in socketids:
        path = '_calc/image_{:d}'.format(socket)
        os.makedirs(path, exist_ok=True)

        with open(path + '/dftb_in.hsd', 'w') as file:
            driver = '  <<+ ../../dftb_in.hsd' + '\n\n' + \
            '+Driver = +Socket ' + '{\n' + \
            '    !File = "dftbplus_{:d}"'.format(socket) + '\n}'
            file.write(driver)


if __name__ == "__main__":
    main()
