
============================================
Error handling with REKS
============================================

************************
Oscillation of energy
************************

The oscillation of SA-REKS energy in SCC process may appear as follows::

    646       -8.9586042136       0.0000229868     0.621642    0.00993712
    647       -8.9586272004      -0.0000229868     0.601626    0.00993712
    648       -8.9586042136       0.0000229868     0.621642    0.00993712
    649       -8.9586272004      -0.0000229868     0.601626    0.00993712
    650       -8.9586042136       0.0000229868     0.621642    0.00993712
    651       -8.9586272004      -0.0000229868     0.601626    0.00993712
    652       -8.9586042136       0.0000229868     0.621642    0.00993712
    653       -8.9586272004      -0.0000229868     0.601626    0.00993712
    654       -8.9586042136       0.0000229868     0.621642    0.00993712

These oscillations mostly arise from the degeneracy of frontier orbitals. As the difference
between frontier orbitals gets smaller, it is hard to determine the order of frontier orbitals
in the SCC process. As a result, oscillating phenomena appear in those cases. To solve this,
the user can increase the shift value by setting ``Shift`` in ``REKS`` block. Then, it may
show correct variational SA-REKS energy in the SCC process::

    111       -8.9586871392      -0.0000000004     0.610476    0.00000117
    112       -8.9586871395      -0.0000000003     0.610474    0.00000110
    113       -8.9586871398      -0.0000000003     0.610473    0.00000103
    114       -8.9586871401      -0.0000000003     0.610471    0.00000097

  --------------------------------------------------
   Final SA-REKS(2,2) energy:      -8.95868714

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -8.97337341   1.220943  0.779057  0.00
    OSS   -8.94400087   1.000000  1.000000  0.00
   Trip   -8.97307659   1.000000  1.000000  1.00
  --------------------------------------------------

Thus, one can check the FONs become :math:`n_a` ~1.2 and :math:`n_b` ~0.8 in this example. As the
orbital energies of frontier orbitals are close each other, the user should increase the shift
value as well as ``MaxSCCIterations``. The oscillation of SA-REKS energy can also occur with
following irregular pattern::

     84      -11.5656112187       0.0000347081     0.535896    0.00170908
     85      -11.5616527044       0.0039585143     0.533173    0.01062304
   fa =  0.471347, MO swap between a(  14) and b(  15) occurs
     86      -11.5390645462       0.0225881582     0.442694    0.05641408
   fa =  0.460378, MO swap between a(  14) and b(  15) occurs
     87      -11.5618504576      -0.0227859114     0.420756    0.00750593
   fa =  0.481839, MO swap between a(  14) and b(  15) occurs
     88      -11.5564810194       0.0053694382     0.463678    0.02941886
     89      -11.4916836033       0.0647974162     0.993664    0.50049487
     90      -11.5138647944      -0.0221811911     0.999927    0.21955679

******************************************
Selection of minimized energy functional
******************************************

A benzene has two degenerate :math:`\pi` orbitals and two degenerate :math:`\pi^*` orbitals. Thus,
if one want to calculate a benzene molecule, then the active space should include two :math:`\pi`
and two :math:`\pi^*` orbitals. Hence, (2,2) active space may be insufficient to investigate this
system, but the ground state can be optimized with DFTB/SSR(2,2). Since there is a ambiguous point
to decide the active space, the OSS state cannot be correctly generated. Thus, only single-state
REKS can be used for a benzene molecule, and the resulting energy of PPS state is given in::

  --------------------------------------------------
    Final REKS(2,2) energy:     -12.56867223

   State     Energy      FON(1)    FON(2)   Spin
    PPS  -12.56867223   2.000000  0.000000  0.00
  --------------------------------------------------

Note that the orbital energy of upper frontier orbital (in this case, 16th orbital) in *band.out* file
is not correctly provided since there is a ambiguous part in single-state REKS method when the Fock
matrix is constructed as follows::

  14    -6.700  2.00000
  15    -6.700  2.00000
  16    -7.583  0.00000
  17    -1.399  0.00000

This problem only occurs in single-state REKS, thus the band energies are correctly calculated for
SA-REKS or SSR method.

According to the system, there may exist a system to show unstability for a lowest excited singlet
state, which corresponds to the OSS state. In this case, including only PPS state to minimized energy
functional is recommended, then the geometry will be optimized successfully. 

***************************
Fail of CP-REKS equations
***************************

In general, preconditioned conjugate-gradient algorithm is used to solve CP-REKS equations in analytic
gradient. It calculates quickly and accurately for most molecules, but it may fail to calculate a
preconditioner due to the symmetry. Its stable geometry of an acetylene is planar so it shows a high
symmetry, thus a singularity appears as follows::

  ERROR!
  -> A singularity exists in preconditioner for PCG, set Preconditioner = No

In this case, conjugate-gradient without preconditioner helps to solve CP-REKS equations. Then, the
CP-REKS equations are now successfully solved as follows::

  ----------------------------------------------------------------------------------
   Solving CP-REKS equation for  1 state vector...
  ----------------------------------------------------------------------------------
    CG solver: Constructing Y initial guess
    CG solver: Iteration    1    eps =    0.000000000000
    Convergence reached in CG solver after    1 iterations
    CG solver: Calculating converged R, Z, Q2 matrix
  ----------------------------------------------------------------------------------

Therefore, the user can check the geometry and symmetry when these errors occur.


