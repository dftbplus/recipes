.. highlight:: none
.. _sec-sockets:

********************
Socket-Communication
********************

[Input: `recipes/interfaces/ase/sockets/`]

For calculations that heavily rely on file-IO, it is better to use communication
between DFTB+ and ASE with the i-PI protocol. In these cases, the reduction of
the wallclock time can be significant.

Note: At this time, only communication with ASE or i-PI has been tested. In
general, however, other software which correctly implements the i-PI protocol
should be supported.

.. note::

    To enable socket-communication in DFTB+ the `WITH_SOCKETS` flag in the
    configuration file `config.cmake` must be set to TRUE, before starting the
    compilation process!

Geometry Optimization by ASE
============================

Providing the input for DFTB+
-----------------------------

To establish a connection via socket communication, DFTB+ is called with the
appropriate ``Socket{}`` driver option::

    Driver = Socket {
      File = "dftbplus"
      Protocol = i-PI {}
      MaxSteps = 1000
      Verbosity = 0
    }

Which instructs DFTB+ to read geometry driving commands via a named temporary
communication file (stored in /tmp/). In this example, due to a bug in the
current ASE implementation of the i-PI protocol (see note below), the code is
asked to perform up to 1000 geometry steps under control of the external
interface then stop. For external driving codes that support the `EXIT` command
in the protocol, instead DFTB+ can continue running until told to stop by
using::

  MaxSteps = -1

Due to the nature of the i-PI protocol, a dummy geometry must be provided, from
which the number and species of the atoms, as well as unit cell vectors will be
read. In general, the specification of arbitrary positions and cell vectors is
possible. However, to avoid errors due to a different ordering of the atoms in
the ASE-geometry and the species in the dummy-geometry, this example uses the
ability of ASE to read in and write out .gen files. The output geometry
(`geo.gen`) is then used as input for DFTB+, avoiding the possibility of the
above mentioned errors::

    Geometry = GenFormat {
      <<< "geo.gen"
    }

In this example, calculation of a water molecule is performed without charge
self-consistency and the required Slater-Koster files are stored in the
top-level directory (same folder as `dftb_in.hsd`)::

    Hamiltonian = DFTB {
      SCC = No

      MaxAngularMomentum = {
        O = "p"
        H = "s"
      }

      SlaterKosterFiles = Type2FileNames {
        Prefix = "./
        Separator = "-"
        Suffix = ".skf"
      }
    }

For the full input file assembled from the code snippets shown here, please
consult the archive, whose location is stated at the begin of this section.

Geometry and main script
------------------------
The geometry used is a simple water molecule, whose file location is specified in the
main script (`GEO_PATH`)::

    3  C
     O H
         1    1    0.00000000000E+00  -0.10000000000E+01   0.00000000000E+00
         2    2    0.00000000000E+00   0.00000000000E+00   0.78306400000E+00
         3    2    0.00000000000E+00   0.00000000000E+00  -0.78306400000E+00

Following the necessary imports, the main script defines the type of socket to
be used (UNIX socket) as well as the path to the DFTB+ executable (`DFTBP_PATH`)
and the geometry (`GEO_PATH`). The main method then first reads in the geometry
and immediately writes it out for the reasons mentioned above. Subsequently,
the ASE driver ``BFGS()`` for geometry optimization is specified and where to
write the trajectory and the logfile.

Finally, the convergent forces and the corresponding energy of the system can
be read out directly by ASE:

.. code-block:: python

    import sys
    from subprocess import Popen
    from ase.io import read, write
    from ase.optimize import BFGS
    from ase.calculators.socketio import SocketIOCalculator


    UNIXSOCKET = 'dftbplus'
    DFTBP_PATH = ''
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

    if __name__ == "__main__":
        main()

.. note::

    To correctly close sockets on the ASE side, call `calc.close()` at the end
    or, more elegantly, enclose the class ``SocketIOCalculator`` using the
    `with` statement as done in the example shown here. Nevertheless, in the
    current state of ASE, the socket gets closed without warning missing the
    'EXIT' string of the i-PI protocol, which always leads to an error message
    issued by DFTB+ at the end of a calculation driven by socket-communication.
