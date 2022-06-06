.. highlight:: none
.. _sec-md-bo:

****************************
Dynamics in the ground state
****************************

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/md/`]

Dynamics on the Born-Oppenheimer ground state energy surface can be performed in
DFTB+ by setting the input geometry driver to be `VelocityVerlet` ::
  
  Driver = VelocityVerlet{
    TimeStep [fs] = 1.0
    Thermostat = NoseHoover {
      Temperature [Kelvin] = 400
      CouplingStrength [cm^-1] = 3200
    }
    Steps = 20000
    MovedAtoms = 1:-1
    MDRestartFrequency = 100
  }

The velocity Verlet driver should have a time step on the scale of ~10x the
highest vibrational period in the system. 1 fs is a common choice. The  

This input file specifies the initial velocities of the atoms
(alternatively they can be generated from a Maxwell-Boltzmann
distribution, see the next section) ::
  
  Velocities [AA/ps] {
    0.63060001     10.71652407      0.41599521
    -4.78167517     -0.67726160      6.81193886
    .
    .
  }

(see :ref:`md-restart` for restarting calculations from a previous MD
simulation)

Thermalising a system
---------------------

.. only :: builder_html or readthedocs

[Input: `recipes/moleculardynamics/thermalise/`]
  
The initial velocities of atoms can be user supplied, however it is more common
to generate them by thermalising the system starting from an initial
Maxwell-Boltzmann distribution of atomic velocities. These can be generated for
example by using the following input::

  Thermostat = NoseHoover {
    # Target temperature
    Temperature [Kelvin] = 400
    # Approximately the highest vibrational frequency of the molecule
    CouplingStrength [cm^-1] = 3200
  }

within the ``VelocityVerlet`` input block. The initial values of the velocities
are set from a random number generator, hence to make the calculation repeatable
this is set to a specific value ::

  Options = {
    RandomSeed = 3871906
  }

However, for real calculations, it would be common to use a fully random choice,
by omitting the ``RandomSeed`` value.

