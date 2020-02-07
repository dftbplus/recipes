.. highlight:: none
.. _sec-interfaces-ase-fileio:

*******
File-IO
*******

For calculations without heavy file-IO, i.e. systems whose wallclock time is
dominated by the electronic structure part of the calculation (due to extensive
SCC cycles or large system size), the communication between DFTB+ and external
software via file-IO may be suitable.

Calling DFTB+ via ASE
=====================

[Input: `recipes/interfaces/ase/fileio/1_conjgrad/`]

In order for ASE to find the DFTB+ executable and Slater-Koster files,
environment variables must be set. Assuming the use of the BASH shell, this is
done as follows (in general an adaptation to your specific computing environment
will be necessary)::

    $ DFTB_PREFIX=~/slakos/mio-0-1/
    $ DFTB_COMMAND=~/dftbplus/bin/dftb+

In this case the geometry optimization of a simple water molecule serves as an
example::

    3  C
     O H
         1    1    0.00000000000E+00  -0.10000000000E+01   0.00000000000E+00
         2    2    0.00000000000E+00   0.00000000000E+00   0.78306400000E+00
         3    2    0.00000000000E+00   0.00000000000E+00  -0.78306400000E+00

Following the necessary imports, the main script defines the path to the DFTB+
executable and the geometry .gen file containing the water molecule shown above.
The main method then first reads in the geometry. Subsequently, a ``DFTB()``
calculator object with options that are crucial for the actual calculation is
then instantiated. Here we perform a self-consistent calculation with the DFTB+
``ConjugateGradient{}`` method as the driver.

The calculator gets attached to the system or geometry respectively and the 
calculation is started.

Finally, the convergent geometry as well as the forces and corresponding energy 
of the system can be read out directly by ASE:

.. code-block:: python

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

    if __name__ == "__main__":
        main()

The script above causes ASE to create an input file (`dftb_in.hsd`) with the
specified options and invokes DFTB+ in the corresponding directory.

Geometry Optimization by ASE
============================

[Input: `recipes/interfaces/ase/fileio/2_bfgs/`]

Apart from the invocation of DFTB+ via file-IO, the use of DFTB+ as an energy/
force engine in conjunction with an external geometry driver is possible. To do 
so, again, set the required environment variables (see explanation above) and 
consider a water molecule of the following geometry::

    3  C
     O H
         1    1    0.00000000000E+00  -0.10000000000E+01   0.00000000000E+00
         2    2    0.00000000000E+00   0.00000000000E+00   0.78306400000E+00
         3    2    0.00000000000E+00   0.00000000000E+00  -0.78306400000E+00

Similar to the previous section, the path to the .gen file containing the 
geometry is defined first, followed by the calculator. Subsequently, the driver 
of the geometry optimization is specified and where to write the trajectory and 
the logfile. In this case, ``BFGS()`` is used, representative of all the 
calculators provided by ASE:

.. code-block:: python

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

    if __name__ == "__main__":
        main()

The script shown causes ASE to generate appropriate input files for each step of
the geometry optimization. Note that this can lead to heavy file-IO and thus a
significant increase in wallclock time, depending on the speed of the storage
used. Additionally, the self-consistency is started from fresh for each
structure, substantially increasing the number of SCC cycles. Therefore it is
advisable to perform such calculations on a ramdisk or even better via
:ref:`sec-sockets`.
