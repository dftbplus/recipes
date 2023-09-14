.. _sec-interfaces-ipi:

****
i-PI
****

In these examples, the `i-PI <http://ipi-code.org/>`_ will be used to
perform an advanced type of molecular dynamics simulation, with DFTB+
suplying the forces and energies, while i-PI externally controls the
simulation.

Originally designed for path-integral molecular dynamics, i-PI [iPI2]_
is a python based program which can communicate with various
electronic structure and interatomic potential codes, using them to
calculate forces and energies under its control. It now includes a
range of advanced dynamics and structural drivers. These include:

* Closed and open path integration to calculate quantum nuclear
  positions and dynamics
* Molecular dynamics in various ensembles (NVE, NVT, ...)
* Thermodynamic integration, replica exchange and umbrella sampling
  to efficiently calculate free energies and sample free energy
  landscapes
* Hessian evaluation for vibrational modes, free energy or instanton
  dynamics
* Geometry optimisation and transition state search methods (nudged
  elastic band, NEB, is temporarily disabled in i-PI 2.0, see
  :ref:`sec-interfaces-ase-neb` for an alternative)

See the full `list of i-PI features
<http://ipi-code.org/about/features/>`_ for more details.

The i-PI code `documentation <https://ipi-code.org/i-pi/index.html>`_
contains install information and comand references. There are examples
in the `repository <https://github.com/i-pi/i-pi>`_ for interfacing
with various codes, and `tutorials
<https://github.com/lab-cosmo/pimd-tutorial>`_.

Here we will use the `conda-forge
<https://anaconda.org/conda-forge/i-pi>`_ packaging.

The i-PI protocol to control DFTB+ is also supported by other driver
codes, such as the :ref:`sec-interfaces-ase`

Input for i-PI
==============

The i-PI code requires a control file writen in XML, along with the
specific input file for the code it is driving. Depending on the
calculation, i-PI can drive multiple instances of a code, hence
allowing for parallel execution. This is particularly important for
cases where there are multiple coupled structures (beads).

i-PI, depending on mode, also requires an xyz formatted file of
coordinates. The structure is assumed to be periodic, with the
boundary conditions typically stored in the comment line of the file,
for example as a cubic box ::

   7
   # CELL(abcABC):   26.45886    26.45886    26.45886    90.00000    90.00000    90.00000
    O      0.57221124     -0.82099183     -0.34085782      6.48095497
    O      2.93098657     -1.00085694      0.03450005      6.48103367
    H      0.11554376     -1.32891000     -1.03000320      0.60545425
    H      1.75146200     -0.91094158     -0.15265030      0.61619753
    H      3.30479374     -1.55123800      0.74081807      0.60545128
    H      3.58125213     -0.53348695     -0.51336512      0.60545835
    H      0.00591684     -0.23123720      0.18168415      0.60544996

The structure can be for a non-periodic system (as far as DFTB+ is
concerned), but i-PI behaves as though it is periodic.

DFTB+ requirements
==================

At compile time, the cmake configuration for DFTB+ must have the the
socket interface enabled::

  -DWITH_SOCKETS=YES

before building the code. The conda DFTB+ images come with this
already enabled.

On starting, each instance of DFTB+ will read the usual `dftb_in.hsd`
file, before connecting to an i-PI server. The `dftb_in.hsd` file
should contain the usual components of the DFTB+ input, including an
initial geometry (which is used by DFTB+ to assign the type of
boundary conditions and chemical species to each atom).

In order to be externally controlled by i-PI, the `Driver {}` for
DFTB+ should be set to receive external commands. This can either be
via a designated file in the `/tmp/` directory of your machine (used
for the examples here) or by connecting to a specified IP address and
port number. The file based communication looks like ::

  Driver = Socket { # communicate via a unix socket file in /tmp
    Protocol = i-PI {} # i-PI interface
    MaxSteps = -1 # run until terminated
    File = "zundel" # name used for the communication file
  }

The method and location of the communication given in `dftb_in.hsd`
should, of course, match the choice made in the i-PI `input.xml` file ::

  <!-- Communicate via the file /tmp/ipi-zundel, do not wrap
  coordinates by periodic boundary conditions -->
  <ffsocket mode='unix' name='zundel' pbc="false">
    <latency> 1.0e-1 </latency>
    <timeout> 6.0e+02 </timeout>
    <address>zundel</address>
  </ffsocket>

Geometry optimisation
=====================

[Input: `recipes/interfaces/ipi/zundel/`]

An example using a Zundel ion (H\ :sub:`5`\ O\ :sub:`2`:sup:`+`) with
on-site corrected DFTB2 embedded in a implicit water solvent is
provided as an example. The implicit solvent lacks the
hydrogen-bonding network that would be present in real water.

Starting a calculation
----------------------

To start the i-pi program, in the directory containing `minimize.xml`
just type ::

  i-pi minimize.xml

Then, potentially in a different directory, start the dftb+
calculation with the `dftb_in.hsd` present at that location ::

  dftb+ > /dev/null

What happens if you start the DFTB+ binary first before i-PI?

The i-PI calculation has very tight tolerances for geometry
optimisation, so will take a few steps, even for this relaxed
structure.

Quantum atomic dynamics
=======================

[Input: `recipes/interfaces/ipi/zundel/`]

Path integral molecular dynamics can be used to sample quantum
behaviour at finite temperatures. It relies on the equivalence between
the thermal ensemble behavior of a set of connected classical systems
and the quantum behaviour of a single system.

If the classical systems are coupled as `ring polymers`, this allows
the determination of equilibrium properties, such as the average
location of atoms and distribution around these positions. The i-PI
code samples the quantum mechanics as distinguishable particles, hence
bosonic or fermionic statistics are not included.

The example file `nvtPI.xml` performs PI-MD at 300 Kelvin using 8
beads and the stochastic velocity-rescaling thermostat (SVR, [svr]_)
::

  i-pi nvtPI.xml

This thermostat is relatively insensitive to parameter choices and
does not strongly affect dynamics.

The temperature profile for the first 400 steps (and the profile if
continued for longer) is shown below

  .. figure:: ../../_figures/interfaces/ipi/zundel/H5O2.png
     :height: 40ex
     :align: center
     :alt: Temperature of simulation up to and after restart.

400 steps is insufficent to thermalise the system to the target
temperature, but the dynamics of the individual beads, and the *path
centroid* (i.e., the average position of the quantum particles) are
output, along with data gathered during the run in the files

  +-----------------------+-----------------+
  | File                  | Contents        |
  +=======================+=================+
  | simulation.xc.xyz     | path centroid   |
  +-----------------------+-----------------+
  | RESTART               | restart data    |
  +-----------------------+-----------------+
  | simulation.out        | collected data  |
  +-----------------------+-----------------+
  | simulation.pos_*.xyz  | individual beads|
  +-----------------------+-----------------+

  
The frequency at which these files are appended is set in `nvtPI.xml`

Stopping programs
-----------------

In the directory where i-PI is running, creating a specific file will
cleanly stop the calculation ::

  touch EXIT

If this file is present, i-PI will halt (or not start in the first
place).

Similarly, in the DFTB+ working directory ::

  touch stop_driver

halts DFTB+. And again, you will need to remove this stop file (if
present) before starting DFTB+.

Likewise, `ctrl + c` will stop either program if issued on the
connected terminal.

Test stopping either the i-PI or DFTB+ binaries, what happens?

Multiple DFTB+ clients
----------------------

A simple bash shell loop to run four separate DFTB+ clients on a
shared file system could look something like ::

  for a in $(seq 4)
  do
    mkdir $a
    cp dftb_in.hsd $a
    cp start.xyz $a
    cd $a
    dftb+ > /dev/null &
    cd ..
  done

(it is probably a good idea to set the shell variable
`OMP_NUM_THREADS=1`).  Each DFTB+ instance is then run in a separate
directory, but communicate with the same i-PI instance via a file in
`/tmp`.

Restarting a calculation
------------------------

[Input: `recipes/interfaces/ipi/zundel/restart_data/`]

To restart i-PI is simple, as it regularly generates checkpoint and
restart files ::

  i-pi RESTART

Then start the DFTB+ client again.

In the case that i-PI completed normally, but you want to continue a
trajectory, you will need to edit the `RESTART` file to increase the
number of steps ::

  <step>40000</step>
  <total_steps>40400</total_steps>

The first line is the step reached when the restart file was written,
while the second line is the required total number of steps.

Continued versions of the output files will then be generated. The
system has been thermalized by this stage, but examine the centroid
structure to see what happened to thew ion and how its motion changes
as it thermalizes.


You can alspo edit a restart file to continue a calculation started as
one type as a new type of evaluation, for example converting an
initial NVE calculation to continue as an `NPT ensemble
<https://ipi-code.org/i-pi/tutorials.html#modifying-the-restart-file>`_
calculation.


More advanced i-PI applications
===============================

Coloured-noise thermostats
--------------------------

To reduce the number of beads required in path-integral simulations,
i-PI can use thermostats that use noise designed to sample high
frequency modes more efficiently.

See the `GLE4MD <http://gle4md.org/>`_ website for more information.

Instantons
----------

For barrier crossings, where quantum effects are important, i-PI also
has the option of applying the instanton semi-classical
approximation. Here, open paths are generated from an initial
structure and hessian. See [iPI2]_ and this `tutorial
<https://github.com/lab-cosmo/pimd-tutorial>`_ for details.

References
==========

.. [iPI2] i-PI 2.0: A universal force engine for advanced molecular
           simulations, V. Kapil et al. Computer Physics Communications (2018)
           DOI: `10.1016/j.cpc.2018.09.020
           <https://doi.org/10.1016/j.cpc.2018.09.020>`_

.. [svr] Canonical sampling through velocity
	 rescaling, G. Bussia, D. Donadio,
	 and M. Parrinello, J. Chem. Phys. (2007) DOI:
	 `10.1063/1.2408420
	 <https://aip.scitation.org/doi/10.1063/1.2408420>`_
