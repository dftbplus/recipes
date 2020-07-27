.. highlight:: none
.. _sec-interfaces-ase-neb-sockets:

****************************
NEB via Socket-Communication
****************************

[Input: `recipes/interfaces/ase/neb/sockets/`]

.. note::

    To enable socket-communication in DFTB+ the `WITH_SOCKETS` flag in the
    configuration file `config.cmake` must be set to TRUE, before starting the
    compilation process!

.. _sec_interfaces_ase_neb_sockets-input:

Providing the input for DFTB+
-----------------------------

The composition of the input is analogous to the explanations in
:ref:`sec-interfaces-ase-sockets`. DFTB+ is invoked with the ``Socket{}``
driver, along the usual specifications for the Slater-Koster files, assuming
that they are stored in the top-level directory of the calculation
(same folder as `dftb_in.hsd`). The name `dftbplus` of the temporary file, which
is required for communication via the i-PI protocol, initially only serves as a
placeholder and is later replaced by the Python script used to run the NEB
calculation. The extended HSD format allows to override properties
(cf. :ref:`sec_interfaces_ase_neb_sockets-mainscript`) that were set earlier and
is well suited for this use case. This effort is being made to ensure that each
image has its own calculator attached::

    Geometry = GenFormat {
        <<< "../../dummy.gen"
    }

    Driver = Socket {
        File = "dftbplus"
        Protocol = i-PI {}
        MaxSteps = -1
        Verbosity = 10
    }

    Hamiltonian = DFTB {
      SCC = Yes
      SCCTolerance = 1.0e-06

      MaxAngularMomentum = {
        H = "s"
        N = "p"
      }

      SlaterKosterFiles = Type2FileNames {
        Prefix = "../../slakos/"
        Separator = "-"
        Suffix = ".skf"
      }
    }

    ParserOptions {
      ParserVersion = 8
    }


.. _sec_interfaces_ase_neb_sockets-mainscript:

Geometry
--------

For the geometries used as the initial and final state, please refer to section
:ref:`NEB (File-IO): Geometry <sec_interfaces_ase_neb_fileio-geometry>`.

Main script
-----------

In many respects, the script at hand is similar to that for File-IO calculations
(cf.
:ref:`NEB (File-IO): Main script <sec_interfaces_ase_neb_fileio-mainscript>`).
Nevertheless, a slight modification is necessary:

Performing the NEB calculation requires setting up a folder structure so that
the calculation of an image is carried out in a separate directory and is
assigned to its own ASE calculator, while communicating via a separate UNIX
socket or temporary file respectively. The Python script shown below is a
suitable and simple implementation to run a calculation with `NIMAGES` as the
number of intermediate images (in this case set to 13).

DFTBP_PATH specifies the path to the DFTB+ executable and HSD_INPUT the path to
the input file discussed above (cf.
:ref:`sec_interfaces_ase_neb_sockets-input`).

A list with the working directories for the DFTB+ calculations is then created
and the input file is read. Method `write_modhsd` creates the separate
subfolders for the different DFTB+ image calculations and copies the input file
into the created directories. Using the capabilities of the extended HSD format,
the filename `dftbplus` of the temporary communication file gets overwritten
according to the corresponding image number.

A list of instantiated socket calculators is created and its elements attached
to the single intermediate images. Sockets and clients are invoked and the
calculators are closed after the calculations have been carried out:

.. code-block:: python

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

Analysis
--------

For a short introduction to the evaluation of the results, please consult the
relevant section
:ref:`NEB (File-IO): Analysis <sec_interfaces_ase_neb_fileio-analysis>`.
