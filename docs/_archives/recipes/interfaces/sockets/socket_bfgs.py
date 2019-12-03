import sys
from subprocess import Popen
from ase.io import read, write
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator


UNIXSOCKET = 'dftbplus'
DFTBP_PATH = 'dftb+'
GEO_PATH = './H2O_cluster.gen'

def main():
    '''Main driver routine.'''

    system = read(GEO_PATH, format='gen')
    write('geo.gen', system, format='gen')

    opt = BFGS(system, trajectory='opt.traj', logfile='opt.log')

    with SocketIOCalculator(log=sys.stdout, unixsocket=UNIXSOCKET) as calc:
        Popen(DFTBP_PATH)
        system.set_calculator(calc)
        opt.run(fmax=1.00E-09)

    forces = system.get_forces()
    energy = system.get_potential_energy()

main()

