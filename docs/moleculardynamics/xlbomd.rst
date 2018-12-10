.. highlight:: none

******************
Speeding up SCC MD
******************

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/XLBOMD/`]

Instead of standard Born-Oppenheimer (BO) molecular dynamics, there are several
alternative methods that can propagate both the atomic positions and electronic
structure simultaneously, often with substantial improvements in computational
speed.

These are either the Car-Parrinello scheme or extended Lagrangian (XL)
Born-Oppenheimer dynamics. DFTB+ supports the second of these methods, which
(when stable) produces SCC equivalent dynamics to conventional Born-Oppenheimer
molecular dynamics but with similar (lower) costs to non-SCC DFTB.

This example uses the fast XL-BOMD method, with a single SCC step for each
geometry once the calculation starts. Here the dynamics of an extended
Lagrangian is used to predict the charges for each time step, but without using
a self-consistent loop. Fast XL-BOMD is potentially less stable, so should only
be used for systems which show good SCC behaviour (molecules without
degeneracies or solids with wide band gaps).

To enable XL dynamics, the ``VelocityVerlet{}`` block should be modified to
contain an extra section::

   XlbomdFast {
       IntegrationSteps = 5
       Scale = 0.5
       TransientSteps = 10
   }

This instructs the code to use 5 time steps to integrate the equations of motion
(EOM) for the system with an initial 10 steps of Born-Oppenheimer dynamics
before starting the XL dynamics. The initial BO steps lead to a more stable
start for the XL dynamics. The scaling factor is used to also increase stability
of the XL-EOMs, and should be in the range of (0..1], where the largest value
which is stable should be used (this can only be determined by testing).

Due to the reduced requirements on SCC during XL-BOMD dynamics, the standard
DFTB forces are usually not sufficiently accurate due to increased
self-consistency errors. Hence an extra command should be added to the
``Hamiltonian {}`` block::

  ForceEvaluation = Dynamics

This enables an additional correction in the forces, which can be used for
``Fermi`` fillings. This correction can also be used in other calculations where
the forces are evaluated with limited self-consistency (for example to speed up
geometry relaxation be relaxing the ``SCCTolerance``).

In the special case of the electronic temperature of the system being set to be
0, there is a faster version of this correction which can be used::

  ForceEvaluation = DynamicsT0
