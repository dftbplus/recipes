.. highlight:: none

****************************
Dynamics in the ground state
****************************

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

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md3.dftb_in.hsd>`

The input file specifies initial velocities of the atoms ::
  
  Velocities [AA/ps] {
    0.63060001     10.71652407      0.41599521
    -4.78167517     -0.67726160      6.81193886
    .
    .
  }
  
The initial velocities can be user suplied, however it is more common to
generate them by thermalising the system starting from an initial
Maxwell-Boltzmann distribution of atomic velocities. These can be generated for
example by using the following input

.. only:: builder_html
   
   See :download:`the full input <../_downloads/md4.dftb_in.hsd>`

