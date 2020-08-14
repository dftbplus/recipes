.. _sec-reks:

============================================
Energy and gradient calculations with REKS
============================================

*******************
Single-state REKS
*******************

[Input: `recipes/reks/single_state_reks/`]

In single-state REKS, only ground state is calculated and it can treat the
state with multireference character. For example, a ethylene molecule
with torsional angle of 90 degrees shows degenerate :math:`{\pi}` and
:math:`{\pi}^*` orbitals, which are correspond to highest occupied molecular
orbital (HOMO) and lowest unoccupied molecular orbital (LUMO), respectively.
In this case the SCC-DFTB calculation shows oscillating phenomena for self-consistent
charges and it cannot obtain a correct variational energy since it cannot consider
the static electronic correlations. In addition, REKS can deal with bond-breaking,
diradicals and magnetic couplings.

To deal with the static electronic correlation, REKS first determines the size for
active orbitals, and only two electrons in two frontier orbitals (2,2) is supported
in current version of `DFTB+ <http://www.dftbplus.org>`_.

.. note:: From the active space, REKS automatically construct possible electronic
   configurations, called microstate (See the paper `JCTC, 2019, 15, 3021-3032.
   <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00132>`_). Total 6 microstates
   are generated in REKS(2,2) calculations, for example, 5th and 6th microstates
   express high-spin configurations in (2,2) space. To correctly describe spin
   configurations except doubly-paired microstates, the ``SpinConstants`` block
   is necessary in ``Hamiltonian`` block. For example, the number of up-electron
   in 5th microstate is larger than that of down-electron. Thus, the atomic spin
   constants are used to describe the spin contributions. Note that the user should
   not include the ``SpinPolarisation`` block, only ``SpinConstants`` block is
   needed in ``Hamiltonian`` block.

With this active space, the single-state REKS can be carried out using the
following input that follows::

  Hamiltonian = DFTB {
    SCC = Yes
    SCCTolerance = 1e-6
    MaxSCCIterations = 10000
    Charge = 0.0
    SpinConstants = {
      ShellResolvedSpin = Yes
      H = { -0.072 }
      C = { -0.031 -0.025 -0.025 -0.023 }
    }
    ...
  }
  REKS = SSR22 {
    Energy = {
      Functional = { "PPS" }
    }
    TargetState = 1
    FonMaxIter = 100
    Shift = 1.0
    VerbosityLevel = 1
  }

Note that the ``REKS = SSR22`` block is a separate block like ``ExcitedStates`` block,
which means it is not included in the ``Hamiltonian`` block. In REKS calculation
``Hamiltonian`` block only controls a detailed Hamiltonian setting such as range-separated
hybrid functional or dispersion corrections. The three singlet states can be generated
from (2,2) active space. These are perfectly spin-paired singlet state (PPS), open-shell
singlet state (OSS) and doubly excited singlet state (DES). Combining these singlet states,
we can determine the singlet states to be included in minimized energy functional. In
single-state REKS case we treats only one singlet state and it becomes PPS state. Thus,
you can only include the PPS state into ``Functional`` block in ``Energy`` block.

The energy for PPS state is expressed by

.. math:: E_{PPS} = \sum_L C_L^{PPS} E_L[\rho_L]

where :math:`C_L^{PPS}` is weighting factors of each microstate for PPS state and :math:`E_L`
is the energy of each microstate. Here the weighting factors are calculated using the
fractional occupation numbers (FONs) of frontier orbitals, thus the energy of PPS state is
minimized with respect to the KS orbitals and FONs in single-state REKS.

With this inputs, the standard output shows information about the energy of PPS state and
the FONs contributing the state. The minimized geometry for ethylene is a planar structure
and the resulting energy and FONs are given in the standard output::

  --------------------------------------------------
    Final REKS(2,2) energy:      -4.89357215

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.89357215   1.999990  0.000010  0.00
  --------------------------------------------------


The energy of PPS state is -4.8936 Hartree and the FONs of frontier orbitals are ~2.0 and ~0.0,
respectively. The orbital character for HOMO is :math:`\pi` and the character for LUMO is
:math:`\pi^*`, so the two electrons are occupied in the HOMO. When the =CH2 begins to rotate into
90 degree, the molecular orbitals change from :math:`\pi`, :math:`\pi^*` orbitals to two localized
single :math:`\pi` orbitals. This means the frontier orbitals are now almost degenerate, and this
is well described with REKS(2,2) calculations. The resulting energy and FONs are given in::

  --------------------------------------------------
    Final REKS(2,2) energy:      -4.77774313

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.77774313   1.000004  0.999996  0.00
  --------------------------------------------------


Now, the energy of PPS state becomes -4.7777 Hartree and the FONs become ~1.0 and ~1.0,
respectively. In addition, the user can check the energy contribution for the PPS state from
*detailed.out* file. The energy for spin contribution is ~0.0 Hartree for planar structure,
while the energy becomes now -0.0188 Hartree for 90 degree rotated structure. Overall, we can
conclude that the planar structure is almost consisted of only first microstate, while the
rotated structure is consisted of several microstates.

*******************
SA-REKS
*******************

[Input: `recipes/reks/sa_reks/`]

Single-state REKS can treat only ground state, thus the vertical excitation energy cannot be
calculated with this method. From the restricted open-shell Kohn-Sham scheme, we can construct
the OSS state which is expressed by

.. math:: E_{OSS} = \sum_L C_L^{OSS} E_L[\rho_L]

where :math:`C_L^{OSS}` is weighting factors of each microstate for OSS state. By minimizing
the energy for the ensemble of PPS and OSS states, we can calculate the vertical excitation
energy between them with state-average REKS (SA-REKS). Again we can calculate the energy of
PPS and OSS states of a ethylene molecule with SA-REKS. The ``REKS`` block has now an additional
``Gradient`` block to calculate the gradient for target state::

  REKS = SSR22 {
    Energy = {
      Functional = { "PPS" "OSS" }
    }
    TargetState = 1
    FonMaxIter = 100
    Shift = 1.0
    Gradient = ConjugateGradient {
      CGmaxIter = 100
      Tolerance = 1.0E-8
      Preconditioner = Yes
      SaveMemory = Yes
    }
    VerbosityLevel = 1
  }

At first the user now have to include the OSS state in ``Functional`` block so that the energy
of OSS state is calculated. The resulting energy and additional information is given by::

  --------------------------------------------------
   Final SA-REKS(2,2) energy:      -4.78921495

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.89357215   1.999990  0.000010  0.00
    OSS   -4.68485776   1.000000  1.000000  0.00
   Trip   -4.73085776   1.000000  1.000000  1.00
  --------------------------------------------------

   Lagrangian Wab:   0.00000000  0.00000000

  --------------------------------------------------
   SSR: 2SI-2SA-REKS(2,2) Hamiltonian matrix
                 PPS           OSS
     PPS    -4.89357215    0.00000000
     OSS     0.00000000   -4.68485776
  --------------------------------------------------

   unrelaxed SA-REKS FONs for S0:  1.999990  0.000010

Here, Final SA-REKS(2,2) energy is the energy of ensemble of PPS and OSS states, which is the quantity
to be minimized in SA-REKS formalism. For the planar structure of a ethylene molecule, the energies of
two states are -4.8936 and -4.6849 Hartree, respectively. The FONs for PPS state are ~2.0 and ~0.0, while
those for OSS state are ~1.0 and ~1.0. In addition, the energy of triplet configuration which corresponds
to 5th or 6th microstate is now given in the standard output. Note that this energy is not the energy
of the triplet state. The user can check the successful convergence by comparing two Lagrangian Wab values.
The two Lagrangian values will become almost same if the energy is converged enough. With the energies and
coupling terms between them, we can construct the 2 :math:`\times` 2 Hamiltonian matrix in the state basis.

After the energy calculation is finished, the gradient for target state, which is PPS state in this example,
is calculated and the final gradient appears at the bottom of the standard output. The keywords in the
``Gradient`` block affect coupled-perturbed REKS (CP-REKS) equations which are used to calculate the gradient
of target state. Here we choose conjugate gradient solver and the keywords ``Preconditioner`` and ``SaveMemory``
are used to accelerate the computational speed of CP-REKS. These two keywords will be switched off depending on
users system.

Similar to the one above, the distorted structure can be calculated using SA-REKS and the results are given
in the following::

  --------------------------------------------------
   Final SA-REKS(2,2) energy:      -4.75910919

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.77753765   1.000000  1.000000  0.00
    OSS   -4.74068073   1.000000  1.000000  0.00
   Trip   -4.77753837   1.000000  1.000000  1.00
  --------------------------------------------------

   Lagrangian Wab:  -0.00004046  0.00004230

  --------------------------------------------------
   SSR: 2SI-2SA-REKS(2,2) Hamiltonian matrix
                 PPS           OSS
     PPS    -4.77753765   -0.00000000
     OSS    -0.00000000   -4.74068073
  --------------------------------------------------

   unrelaxed SSR FONs for S0:  1.000000  1.000000

Now the energy of OSS state is -4.7407 Hartree and the coupling between PPS and OSS states becomes non-zero value.
The energy of triplet configuration is similar with the energy of PPS state. Since REKS can consider the static
electronic correlations, it can show a correct shape for the potential energy curve with respect to the torsional
angle of C=C bond. If you want to calculate the energy of DES state, then ``IncludeAllStates = Yes`` keyword in
the ``Energy`` block will show the energy of DES state as well as the FONs.

*******************
SI-SA-REKS
*******************

[Input: `recipes/reks/si_sa_reks/`]

State-interaction SA-REKS (SI-SA-REKS, briefly SSR) energies are obtained by solveing 2 :math:`\times` 2 secular
equation with the possible couplings between the electronic states.

.. math:: \left(\begin{array}{cc} E^{PPS} & \Delta^{SA} \\ \Delta^{SA} & E^{OSS} \end{array}\right)
          \left(\begin{array}{cc} a_{00} & a_{01} \\ a_{10} & a_{11} \end{array}\right) =
          \left(\begin{array}{cc} E^{SSR}_0 & 0 \\ 0 & E^{SSR}_1 \end{array}\right)
          \left(\begin{array}{cc} a_{00} & a_{01} \\ a_{10} & a_{11} \end{array}\right)

By considering the state-interaction terms, SSR states become more reliable states when the excited states are
included. SSR states can be calculated with the ``StateInteractions = Yes`` in ``Energy`` block. For the
planar structure, the resulting energies are given by::

  ----------------------------------------------------------------
   SSR: 2SI-2SA-REKS(2,2) states
                      E_n       C_{PPS}    C_{OSS}
   SSR state  1   -4.89357215  -1.000000   0.000000
   SSR state  2   -4.68485776  -0.000000  -1.000000
  ----------------------------------------------------------------

In this case the ground state is consisted of PPS state, while the lowest excited state is consisted of
OSS state. As the coupling term increases, the difference between SA-REKS and SSR energies becomes larger.


