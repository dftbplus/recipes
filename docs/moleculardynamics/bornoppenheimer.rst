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
   
   See :download:`the full input <../inputs/moleculardynamics/bomd/dftb_in.hsd>`
   and :download:`geometry <../inputs/moleculardynamics/bomd/geom.gen>`.
