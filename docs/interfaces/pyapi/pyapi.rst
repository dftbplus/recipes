.. highlight:: none
.. _sec-interfaces-pyapi:

****************
Python Interface
****************

The Python interface enables DFTB calculations to be performed by external
programs. Results such as system energy, atomic forces and Mulliken charges can
be extracted directly, along with performing tasks like defining calculations in
external potentials or updating the geometry. With the pythonapi (version 0.1)
an existing input file (`dftb_in.hsd`), which contains settings in addition to
the geometry, is required for the initialization of DFTB+.

Regarding the **units**: Just like DFTB+, the interface expects and delivers
values in **atomic units**!

Setting up DFTB+
================

For this special use case, DFTB+ needs to be compiled as a shared library with
API support enabled. At this point, a basic understanding of how to build DFTB+
is assumed (all necessary steps are explained in the INSTALL.rst file in the top
level directory of the DFTB+ repository). The CMake configuration for this case
can be done by settining the WITH_API and BUILD_SHARED_LIBS flags to be TRUE in
the file config.cmake, before starting the configuration and compilation
process. Alternatively, if you do not want to modify files, a construct like the
following is a convenient way to specify these flags on the command line while
configuring with CMake:

``cmake -DBUILD_SHARED_LIBS=1 -DWITH_API=1 ..``

If only the pure compilation process is carried out, the resulting library is
located inside the buid directory at `prog/dftb+/`. If ``make install`` was
executed in addition, there is also a copy placed in the `CMAKE_INSTALL_PREFIX`
path, which defaults to `install/lib/` inside the build directory. Within these
locations, the library files `libdftbplus.*` will be installed.

The path to the resulting shared library must be passed to the interface, so you
will need to know where the `libdftbplus.so` file ends up.

Setting up the interface
========================

You can install the script package via the standard 'python setup' mechanism. If
you want to install it system-wide into your normal python installation, with
the appropriate permissions you can simply issue ``python setup.py`` in the
directory `tools/pythonapi/`. Alternatively, to install it locally in your home
space, use ``python setup.py install --user``. If the local Python installation
directory is not in your PATH, you should add it accordingly.

.. _sec-interfaces-pyapi-input:

Providing the input for DFTB+
=============================

[Input: `recipes/interfaces/pyapi/`]

To initialize an instance of the dftbplus calculator, you will need an initial
file formatted according to the HSD structure for the usual `dftb_in.hsd` input
of DFTB+.

Although the geometry of the structure to be calculated can be later replaced by
your main script, it is necessary first to supply a dummy geometry (suitable for
the desired number of atoms and boundary conditions). It is important to ensure
that the total number of coordinates matches the number of atoms of the real
geometry and that any specified lattice vectors are linearly independent,
otherwise DFTB+ will return an error message. In the case of the supplied |TiO2|
example the following geometry block of `dftb_in.hsd` would be appropriate::

    Geometry = GenFormat {

	6  S
     Ti  O
	1 1    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	2 1    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	3 2    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	4 2    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	5 2    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	6 2    0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
	0.0000000000E+00    0.0000000000E+00    0.0000000000E+00
	1.0000000000E+02    0.0000000000E+00    0.0000000000E+00
	0.0000000000E+00    1.0000000000E+02    0.0000000000E+00
	0.0000000000E+00    0.0000000000E+00    1.0000000000E+02

    }

The section ``Hamiltonian`` contains the usual entries for the calculation
of a |TiO2| supercell::

    Hamiltonian = DFTB {

      Scc = Yes
      MaxSccIterations = 100
      SccTolerance = 1.00e-05

      SlaterKosterFiles = Type2FileNames {
	Prefix = "./slakos/mio-ext/"
	Separator = "-"
	Suffix = ".skf"
      }

      MaxAngularMomentum {
	Ti = "d"
	O  = "p"
      }

      KPointsAndWeights = SupercellFolding {
	4   0   0
	0   4   0
	0   0   4
	0.5 0.5 0.5
      }

    }

To be able to also extract the forces via the interface, it is advisable to
specify the corresponding entry in the ``Analysis{}`` block. Otherwise DFTB+
will issue an error message, asking the user to do so. To ensure backwards
compatibility of the input, the parser version should also be specified::

    Analysis {
      CalculateForces = Yes
    }

    ParserOptions {
      ParserVersion = 8
    }


.. |TiO2| replace:: TiO\ :sub:`2`\

.. _sec-interfaces-pyapi-mainscript:

Main script
===========

The script shown here serves to illustrate the use of the Python interface,
based on the calculation of |TiO2|.

In order to be able to use the interface, the package `dftbplus` must be
imported as the first step. The path `LIB_PATH` to the DFTB+ shared library is
defined (note that the name prefix of the libray file name should also be part
of the path), as well as conversion factors to convert the atom coordinates we
will list in from Ångström into the atomic units (Bohr) required by the
interface.

.. code-block:: python

    import numpy as np
    import dftbplus


    LIB_PATH = '/home/user/libdftbplus'

    # DFTB+ conversion factors
    # (according to prog/dftb+/lib_common/constants.F90)
    BOHR__AA = 0.529177249
    AA__BOHR = 1 / BOHR__AA

At the beginning of the ``main()`` function, the atom coordinates and lattice
vectors are defined. In this case, a conversion to atomic units is necessary,
since a `.gen` block is used whose values are usually in units of Ångström.

.. _sec-interfaces-pyapi-codeblock1:

.. code-block:: python

    def main():
	'''Main driver routine.'''

	# coordinates of TiO2, in Ångström
	coords = np.array([
	    [-0.016726922839251,  0.016725329441158, -0.000003204152532],
	    [-0.016726505918979,  1.920201169305565, -7.297102897292027],
	    [ 0.017412997824265, -0.024318617967798,  2.005339137853385],
	    [ 1.920770753428742, -0.024319922392223, -4.437737763954652],
	    [ 0.024319174400169, -0.017404302527510, -2.005347277168561],
	    [ 0.024317270342179,  1.886164739806594, -5.291732430733527]])

	# lattice vectors of TiO2, in Ångström
	latvecs = np.array([
	    [-1.903471721000000,  1.903471721000000,  4.864738245000000],
	    [ 1.903471721000000, -1.903471721000000,  4.864738245000000],
	    [ 1.903471721000000,  1.903471721000000, -4.864738245000000]])

	# conversion to atomic units
	coords *= AA__BOHR
	latvecs *= AA__BOHR

An object of the DftbPlus class is instantiated, which requires the location of
the shared library `libpath`, the HSD input file `hsdpath` and the name of the
log file `logfile` to be optionally specified. These keywords have default
values './libdftbplus', './dftb_in.hsd' and None if not set explicitly.  Note,
that adding the shared library extension to `libpath` is not essential.  Since
the extension can be system dependent, it is guessed by the interface if
missing. If logfile=None is specified, the output of the calculation gets
printed to stdout.

After instantiation, the geometry can set or replaced; for periodic structures,
lattice vectors can be specified in addition to the absolute coordinates.

The DFTB+ calculations are carried out automatically, as soon as the
corresponding get_* methods are called. To correctly finalize the DFTB+ object,
use the ``close()`` method.

.. _sec-interfaces-pyapi-codeblock2:

.. code-block:: python

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
	grosschg = cdftb.get_gross_charges()

	# finalize DFTB+ and clean up
	cdftb.close()


    if __name__ == "__main__":
	main()

As always, please consult the archive to obtain the complete, connected script.
To do so, follow the path mentioned above.

The Interface is also capable of defining a population (in)dependent external
potential. This is covered in the following two sections
(:ref:`extpot <sec-interfaces-pyapi-extpot>`,
:ref:`qdepextpot <sec-interfaces-pyapi-qdepextpot>`).

.. _sec-interfaces-pyapi-extpot:

Using a population independent external potential
=================================================

[Input: `recipes/interfaces/pyapi/extpot/`]

A external potential which does not depend on the Mulliken charges in the
calculation can be included with only a small addition in the script. The DFTB+
object has a method ``set_external_potential()``, which should be relatively
self-explanatory. The external potential at the position of the QM-atoms is
given as a positional argument. If forces are required, the gradient of the
external potential at each atom can be passed additionally as the keyword
argument `extpotgrad`.

Therefore, after the initialization of the DFTB+ object, the following code is
inserted:

.. code-block:: python

    # example values of extpot and extpotgrad used here were
    # taken from file: test/api/mm/testers/test_extpot.f90
    extpot = np.array([-0.025850198503435,
                       -0.005996294763958,
                       -0.022919371690684])

    extpotgrad = np.array([
        [0.035702717378527,  0.011677956375860, 0.009766745155626],
        [0.023243271928971, -0.000046945156575, 0.004850533043745],
        [0.016384005706180,  0.004608295375551, 0.005401080774962]])

    # set external potential and its gradients
    cdftb.set_external_potential(extpot, extpotgrad=extpotgrad)

.. _sec-interfaces-pyapi-qdepextpot:

Population dependent external potential
=======================================

[Input: `recipes/interfaces/pyapi/qdepextpot/`]

This section deals with the capability of the interface to run calculations with
a population dependent external potential, i.e. arrising in cases like
polarizable surroundings where the applied field responds to the state of the QM
calculation. Since in general only the user knows how to calculate this type of
potential, callback functions can be defined which will then be executed at runtime.

The DFTB+ object provides a method ``register_ext_pot_generator()`` that takes
care of the registration of the callback functions. As the first positional
argument of this method, an arbitrary pointer can be specified. DFTB+ will pass
back this pointer unaltered when calling the registered functions. You can
typically use it to pass a pointer to the data or a Python object (class) which
contains the necessary data for the potential calculation. If your data is in
the global space and you do not need it, pass None (or equivalent). The second
and third positional arguments have to be the function that provides the
external potential and its gradients.

Furthermore, the auxiliary class ``PotentialCalculator`` is defined to perform
the actual calculation of the external potential and its gradients. The
structure of a script required for the calculation is explained below, using a
trivial example in which the external potential and gradient are assumed to be
zero. Therefore, this should not change the results of the calculation.

.. code-block:: python

    import numpy as np
    import dftbplus


    LIB_PATH = '/home/user/libdftbplus'


    class PotentialCalculator:
	'''

	   Auxiliary class for calculating the population dependent external
	   potential and its gradients. An instance of this class gets handed over
	   to DFTB+ via the ctypes interface, to handle the necessary callbacks.

	'''


	def __init__(self, qmcoords, mmcoords, mmcharges):
	    '''Initializes a PotentialCalculator object.

	    Args:

		qmcoords (2darray): coordinates of QM-atoms
		    (shape: [qmatoms, 3])
		mmcoords (2darray): coordinates of MM-atoms
		    (shape: [mmatoms, 3])
		mmcharges (1darray): charges of MM-atoms
		    (shape: [mmatoms, 1])

	    '''

	    self._qmcoords = qmcoords
	    self._mmcoords = mmcoords

	    self._qmatoms = np.shape(self._qmcoords)[0]
	    self._mmatoms = np.shape(self._mmcoords)[0]

	    self._mmcharges = mmcharges


	def calc_extpot(self, dqatom):
	    '''Calculates the current external potential
	       using the properties of the MM- and QM-atoms.

	    Args:

		dqatom (1darray): population difference with respect to
		    reference population (usually the neutral atom)
		    Note: population means electrons, so a
		    positive number indicates electron excess

	    Returns:

		extpot (1darray): updated external potential
		    at the position of each QM-atom

	    '''

	    # Note: Some types of potential require knowledge of the
	    # current atomic populations, which is provided by dqatom.

	    extpot = np.zeros(self._qmatoms)

	    return extpot


	def calc_extpotgrad(self, dqatom):
	    '''Calculates the current gradients of the external
	       potential using the properties of the MM- and QM-atoms.

	    Args:

		dqatom (1darray): population difference with respect to
		    reference population (usually the neutral atom)
		    Note: population means electrons, so a ositive number
		    indicates electron excess

	    Returns:

		extpotgrad (2darray): updated potential gradient
		    at the position of each QM-atom

	    '''

	    # Note: Some types of potential require knowledge of the
	    # current atomic populations, which is provided by dqatom.

	    extpotgrad = np.zeros((self._qmatoms, 3))

	    return extpotgrad


    def get_extpot(potcalc, dqatom, extpotatom):
	'''Queries the external potential.

	Args:

	    potcalc (pyobject): instance of a class that provides methods for
	        calculating the external potential and its gradients
	    dqatom (1darray): population difference with respect to reference
	        population (usually the neutral atom)
		Note: population means electrons, so a positive number indicates
		electron excess
	    extpotatom (1darray): potential at the position of each QM-atom
	        Note: it should be the potential as felt by an electron
		(negative potential value means attraction for an electron)

	'''

	extpotatom[:] = potcalc.calc_extpot(dqatom)


    def get_extpotgrad(potcalc, dqatom, extpotatomgrad):
	'''Queries the external potentials gradients.

	Args:

	    potcalc (pyobject): instance of a class that provides methods for
	        calculating the external potential and its gradients
	    dqatom (1darray): population difference with respect to referenc
	        population (usually the neutral atom)
		Note: population means electrons, so a positive number indicates
		electron excess
	    extpotatomgrad (2darray): potential gradient at the position of each
	        QM-atom
		Note: it should be the gradient of the potential as felt by an
		electron (negative potential value means attraction for an
		electron)
	'''

	extpotatomgrad[:, :] = potcalc.calc_extpotgrad(dqatom)

The initialization of the calculator and the definition of the geometry is
completely analogous to the
:ref:`above explanations <sec-interfaces-pyapi-codeblock2>`. Only the
registration of the callback functions is still missing:

.. code-block:: python

    # register callback functions for a qdepextpot calculation
    cdftb.register_ext_pot_generator(potcalc, get_extpot, get_extpotgrad)

Please consult the associated archive with this tutorial to obtain the full
corresponding example.
