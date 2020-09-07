.. _sec-reks:

******************************************
Energy and gradient calculations with REKS
******************************************

.. _single_state_REKS:

Single-state REKS
=================

[Input: `recipes/reks/single_state_reks/`]

In single-state REKS, the ground state is calculated including multi-reference
character. For example, a ethylene molecule with a torsional angle of 90 degrees
has degenerate :math:`{\pi}` and :math:`{\pi}^*` orbitals, corresponding to the
highest occupied molecular orbital (HOMO) and lowest unoccupied molecular
orbital (LUMO) respectively.  In this case a standard SCC-DFTB calculation can
show oscillating phenomena when trying to obtain self-consistent charges and is
missing the static electronic correlations required for a correct variational
energy. Similarly, REKS can deal with bond-breaking, di-radicals and magnetic
couplings which are difficult in standard calculations.

To deal with the static electronic correlation, REKS first determines the size
of the space of active orbitals. Only two electrons in two frontier orbitals,
`(2,2)`, is supported in current version of `DFTB+ <http://www.dftbplus.org>`_.

.. note:: From the active space, REKS automatically construct possible
   electronic configurations, called microstate (See the paper `JCTC, 2019, 15,
   3021-3032.  <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00132>`_). A total
   of 6 microstates are generated in a REKS(2,2) calculations, for example the
   5\ :sup:`th` and 6\ :sup:`th` microstates express high-spin configurations in
   the (2,2) space. To correctly describe the spin polarised configurations in
   the space the ``SpinConstants`` block is required in the ``Hamiltonian``. For
   example, there are more up-electrons in the 5\ :sup:`th` microstate than
   down-electrons, thus atomic spin constants are needed to describe spin
   contributions. Note that the user should not include the ``SpinPolarisation``
   block, only the ``SpinConstants`` block is needed.

With this active space, single-state REKS can be carried out using the following
input::

  Hamiltonian = DFTB {
    SCC = Yes
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

Note that the ``REKS = SSR22`` block is a separate environment in the
calculation, so is not included in the ``Hamiltonian`` block. In a REKS
calculation, the ``Hamiltonian`` only controls detailed Hamiltonian settings
such as use of a range-separated hybrid functional or dispersion
corrections. Three singlet states can be generated from a (2,2) active
space. These are the perfectly spin-paired singlet state (PPS), an open-shell
singlet state (OSS) and the doubly excited singlet state (DES). Combining these
singlet states, we can find the many-particle singlet state which minimises the
energy functional. In the single-state REKS case we treat only one singlet state
of PPS symmetry, thus this is the only possibility that can be selected in the
``Functional`` block of the ``Energy`` in this case.

The energy for the PPS state is expressed as

.. math:: E_{PPS} = \sum_L C_L^{PPS} E_L[\rho_L]

where :math:`C_L^{PPS}` are the weighting factors of each microstate in the
resulting PPS state and :math:`E_L` are their energies. Here the weighting
factors are calculated using the fractional occupation numbers (FONs) of
frontier orbitals, thus the energy of the resulting PPS state is minimised with
respect to both the Kohn-Sham orbitals and the FONs in single-state REKS.

[Input: `recipes/reks/single_state_reks/dftb_in.hsd-0deg`]

With this inputs, the standard output shows information about the energy of the
PPS state and the FONs contributing to it. The minimised geometry for ethylene
is a planar structure and the resulting energy and FONs are given in the
standard output::

  --------------------------------------------------
    Final REKS(2,2) energy:      -4.89357215

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.89357215   1.999990  0.000010  0.00
  --------------------------------------------------

The energy of PPS state is -4.8936 Hartree and the FONs of frontier orbitals are
~2.0 and ~0.0, respectively. The orbital character for HOMO is :math:`\pi` and
the character for LUMO is :math:`\pi^*`, so the two electrons are occupied in
the HOMO.

[Input: `recipes/reks/single_state_reks/dftb_in.hsd-90deg`]

When the torsional angle begins to rotate to 90 degree, the molecular orbitals
change from :math:`\pi` and :math:`\pi^*` orbitals into two localised single
:math:`\pi` orbitals. This means the frontier orbitals are now almost
degenerate, and this is well described with REKS(2,2) calculations. The
resulting energy and FONs are::

  --------------------------------------------------
    Final REKS(2,2) energy:      -4.77774313

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.77774313   1.000004  0.999996  0.00
  --------------------------------------------------

Now, the energy of PPS state becomes -4.7777 Hartree and the FONs become ~1.0
and ~1.0, respectively. In addition, the user can check the energy contributions
to the PPS state from *detailed.out* file. The energy from the spin contribution
is ~0.0 Hartree for the planar structure, while the energy becomes -0.0188
Hartree for the 90 degree rotated structure. Overall, we can conclude that the
planar structure almost consists of only the first microstate, while the rotated
structure consists of several microstates.

.. _sa_reks:

SA-REKS
=======

[Input: `recipes/reks/sa_reks/`]

Single-state REKS can treat only the ground state, thus the vertical excitation
energy cannot be calculated with this method. From the restricted open-shell
Kohn-Sham scheme, we can construct the OSS state, which is expressed as

.. math:: E_{OSS} = \sum_L C_L^{OSS} E_L[\rho_L]

where :math:`C_L^{OSS}` is weighting factors of each microstate making up the
OSS state. By minimising the energy for the ensemble of PPS and OSS states, we
can calculate the vertical excitation energy between them with state-average
REKS (SA-REKS). Again we calculate the energy of the PPS and OSS states of the
two geometries of an ethylene molecule with SA-REKS. The ``REKS`` block has now
an additional ``Gradient`` block to calculate the gradient and optimise for the
target state::

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

The user also now has to include the OSS state in addition to the PPS in the
``Functional`` block so that the energy of both are now calculated. The
resulting energy and additional information is given by::

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

The final SA-REKS(2,2) energy is from the ensemble of the PPS and OSS states,
the energy of which is being minimised. For the planar structure of an ethylene
molecule, the energies of two states are -4.8936 and -4.6849 Hartree,
respectively. The FONs for the PPS state are ~2.0 and ~0.0, while those for OSS
state are ~1.0 and ~1.0. In addition, the energy of the triplet configuration
which corresponds to contributions from the 5\ :sup:`th` and 6\ :sup:`th`
microstates is now given in the standard output. Note that this energy is *not*
the energy of the triplet state. The user can check for successful convergence
by comparing the two Lagrangian *Wab* values, these will become almost the same
if the energy is well converged.

After the energy calculation is finished, the gradient for the target state
(``TargetState = 1``, which is the PPS state in this example) is calculated and the
final gradient appears at the bottom of the standard output. The keywords in the
``Gradient`` block affect coupled-perturbed REKS (CP-REKS) equations which are
used to calculate the gradient of target state. Here we choose a conjugate
gradient solver and the keywords ``Preconditioner`` and ``SaveMemory`` are used
to accelerate the computational speed of CP-REKS. These two keywords may be
switched off depending on the user's system.

Similar to the case above, the distorted structure can be calculated using
SA-REKS and the results are given in the following::

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

Now the energy of OSS state is -4.7407 Hartree. The energy of the triplet
configuration is similar to the energy of PPS state. Since REKS can consider the
static electronic correlations, it can produce the correct shape for the
potential energy curve with respect to the torsional angle of C=C bond. If you
want to calculate the energy of the DES state, then ``IncludeAllStates = Yes``
keyword in the ``Energy`` block will produce the energy of DES state as well as
its FONs.

.. _si_sa_reks:

SI-SA-REKS
==========

[Input: `recipes/reks/si_sa_reks/`]

State-interaction SA-REKS (SI-SA-REKS, briefly `SSR`) energies are obtained by
solving a 2 :math:`\times` 2 secular equation with the possible couplings
between the SA-REKS states.

.. math:: \left(\begin{array}{cc} E^{PPS} & \Delta^{SA} \\ \Delta^{SA} & E^{OSS} \end{array}\right)
          \left(\begin{array}{cc} a_{00} & a_{01} \\ a_{10} & a_{11} \end{array}\right) =
          \left(\begin{array}{cc} E^{SSR}_0 & 0 \\ 0 & E^{SSR}_1 \end{array}\right)
          \left(\begin{array}{cc} a_{00} & a_{01} \\ a_{10} & a_{11} \end{array}\right)

By considering the state-interaction terms, the SSR states become more reliable
when the excited states are included. The SSR states can be calculated with the
``StateInteractions = Yes`` in ``Energy`` block. For the planar structure, the
resulting energies are given by::

  ----------------------------------------------------------------
   SSR: 2SI-2SA-REKS(2,2) states
                      E_n       C_{PPS}    C_{OSS}
   SSR state  1   -4.89357215  -1.000000   0.000000
   SSR state  2   -4.68485776  -0.000000  -1.000000
  ----------------------------------------------------------------

In this case, the ground state consists of the PPS state, while the lowest
excited state is of OSS type. As the coupling term increases, the difference
between the SA-REKS and SSR energies becomes larger.
