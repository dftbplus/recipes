.. _sec-advanced_REKS:

*********************
Advanced calculations
*********************


Non-adiabatic couplings
=======================

[Input: `recipes/reks/nonadiabatic_coupling/`]

The non-adiabatic coupling vectors between the states can be calculated with the
SSR method. By setting ``NonAdiabaticCoupling = Yes`` in ``REKS`` block, the
user can easily calculate the non-adiabatic coupling vectors. When this option is
turned on, SSR shows the gradients for all states and possible coupling
vectors::

  --------------------------------------------------
   Gradient Information
  --------------------------------------------------
   1 st state (SSR)
       0.09734680    -0.00000000     0.00000000
      -0.09734682    -0.00000000    -0.00000000
       0.00639517     0.00417744    -0.00000000
       0.00639517    -0.00417744     0.00000000
      -0.00639516     0.00417743    -0.00000000
      -0.00639516    -0.00417743     0.00000000

   2 st state (SSR)
      -0.12803706    -0.00000000    -0.00000000
       0.12803704    -0.00000000     0.00000000
       0.00639486     0.00417837    -0.00000000
       0.00639486    -0.00417837     0.00000000
      -0.00639485     0.00417836    -0.00000000
      -0.00639485    -0.00417836     0.00000000
  --------------------------------------------------
   Coupling Information
  --------------------------------------------------
   between  1 and  2 states
  --------------------------------------------------
   g vector - difference gradient
       0.11269193    -0.00000000     0.00000000
      -0.11269193     0.00000000    -0.00000000
       0.00000016    -0.00000046     0.00000000
       0.00000016     0.00000046     0.00000000
      -0.00000016    -0.00000046     0.00000000
      -0.00000016     0.00000046     0.00000000

   h vector - derivative coupling
      -0.00031874    -0.00000000     0.00000000
      -0.00031874     0.00000000    -0.00000000
       0.00015937     0.00007242     0.00000000
       0.00015937    -0.00007242     0.00000000
       0.00015937    -0.00007242     0.00000000
       0.00015937     0.00007242     0.00000000

   G vector - GDV
       0.11269193    -0.00000000     0.00000000
      -0.11269193     0.00000000    -0.00000000
       0.00000016    -0.00000046     0.00000000
       0.00000016     0.00000046     0.00000000
      -0.00000016    -0.00000046     0.00000000
      -0.00000016     0.00000046     0.00000000

   H vector - DCV - non-adiabatic coupling
      -0.00152718    -0.00000000     0.00000000
      -0.00152718     0.00000000    -0.00000000
       0.00076359     0.00034700     0.00000000
       0.00076359    -0.00034700     0.00000000
       0.00076359    -0.00034700     0.00000000
       0.00076359     0.00034700     0.00000000
  --------------------------------------------------

The above gradients and coupling vectors are obtained with the planar structure
of a ethylene molecule.  The *g* vector is defined as gradient difference
vector, thus it can be calculated from the difference of SA-REKS
gradients. Similarly to this, *G* vector is calculated from the difference of
SSR gradients. The *h* vector is defined as coupling gradient, so it can be
simply calculated from the gradient of state-interaction terms. The *H* vector
is defined as the derivative of the coupling vectors, thus its norm increases as
the energy gap becomes smaller.

The *g* and *h* vectors can be regarded as the vectors defined through diabatic
states, and the *G* and *H* vectors are defined through the adiabatic (SSR)
states. In general, the non-adiabatic coupling vectors can be used for surface
hopping molecular dynamics, while the *g* and *h* vectors can be used for
minimum energy conical intersection (MECI) optimisation.

Relaxed Density
===============

[Input: `recipes/reks/relaxed_density/`]

The relaxed density is calculated for the ``TargetState`` when the user sets
``RelaxedDensity`` to ``Yes``. The calculation of a relaxed density requires the
information about gradient, thus it can be calculated when the input enables
gradient calculation. When this option is turned on, the relaxed FONs are given
in the bottom of standard output and the *relaxed_charge.dat* file is
generated. It includes the total charge as well as the Mulliken charges of each
atom for target state.  This example shows relaxed charges for the ground
state::

  total charge:     -0.00000000 (e)

  relaxed S0 atomic charge (e)
     atom        charge
        1      -0.18168616
        2      -0.18168616
        3       0.09084308
        4       0.09084308
        5       0.09084308
        6       0.09084308

If one want to exploit SSR with a QM/MM approach, the ``RelaxedDensity`` option
should be used to obtain the gradient of external point charges. Then, the
gradient of external point charges for target state are given in *detailed.out*
file. Therefore, this option can be used to run surface hopping dynamics with a
QM/MM approach.

Spin tuning constants
=====================

The DFTB/SSR method well describes equilibrium geometries and vertical
excitation energies as compared with SSR/wPBEh results. (See the paper `JCTC,
2019, 15, 3021-3032.  <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00132>`_)
However, the behaviours at MECI points sometimes does not match those obtained
with SSR/wPBEh. For example, the n/:math:`\pi^*` type MECI geometry of ethylene
or methyliminium molecule cannot be located with DFTB/SSR, with an incorrect
description of the relative stability of the PPS and OSS states being mostly
responsible. Their relative stability depends on the splitting between the
open-shell singlet microstates and the triplet microstates in the PPS and OSS
energies.

In principle, DFTB/SSR employs spin-polarised DFTB formalism, in which the
spin-polarisation contribution is obtained from the second-order expansion of
the magnetisation density with respect to zero magnetisation electronic
structure. At the n/:math:`\pi^*` type MECI geometries, both frontier orbitals
are located on the same atom. In such a case, the second-order expansion of
magnetisation may not be suitable for the triplet microstates, as the spin
density becomes too large. As a simple solution, the stability between the PPS
and OSS states can be adjusted by scaling the atomic spin constants. For most
molecules the FONs for the PPS state become :math:`n_a` ~ 2.0 and :math:`n_b` ~
0.0, hence the energy of the PPS state is determined by the 1\ :sup:`st`
microstate alone and it is only the energy of the OSS state that depends on the
atomic spin constants.

If the user runs the test calculation included with the main DFTB+ repository in
the `test/prog/dftb+/reks/PSB3_2SSR_rangesep_tuning` directory, the following
results are given in the standard output::

  ----------------------------------------------------------------
   SSR: 2SI-2SA-REKS(2,2) states
                      E_n       C_{PPS}    C_{OSS}
   SSR state  1  -16.39950035  -0.955182  -0.296020
   SSR state  2  -16.38979921   0.296020  -0.955182
  ----------------------------------------------------------------

  H vector - DCV - non-adiabatic coupling
     -0.11443527     0.16066530     0.18737790
      0.15477051    -0.03677914     0.02047645
     -0.81366262     0.56406749    -1.41489688
      0.69877935    -2.68846832     1.70292281
     -0.34893713     0.81739627     0.50771020
     -0.08309257     0.13356401     0.02704612
      0.99884931     0.52013559    -1.04847579
     -0.46122482     0.81158082    -0.36230963
     -0.03763611     0.00213152     0.02451782
      0.65298995     0.42710364    -0.20222417
     -0.52172598    -0.75236687     1.00033386
     -0.24407958    -0.11018546    -0.50337767
      0.04971872    -0.03074284     0.01547016
      0.06968625     0.18189801     0.04542884

It shows the energies at MECI point of a PSB3 molecule, thus the non-adiabatic
coupling vectors show large elements near to the centre C=C bond. The atomic
spin constants can be modified by using ``SpinTuning`` keyword in ``REKS`` block
as follows::

  Reks = SSR22 {
    Energy = {
      Functional = { "PPS" "OSS" }
      StateInteractions = Yes
    }
    TargetState = 2
    FonMaxIter = 30
    shift = 0.3
    SpinTuning = { 3.2 3.2 3.2 }
    TransitionDipole = Yes
    Gradient = ConjugateGradient {
      CGmaxIter = 100
      Tolerance = 1.0E-8
      Preconditioner = Yes
      SaveMemory = Yes
    }
    RelaxedDensity = Yes
    NonAdiabaticCoupling = Yes
    VerbosityLevel = 1
  }

Microstate calculation
======================

[Input: `recipes/reks/microstate/`]

Obviously, the SSR method treats only singlet states like PPS or OSS. If one
want to compare the energy of singlet and triplet states, SSR provides the
energy of a triplet configuration as an alternative which corresponds to the 5\
:sup:`th` or 6\ :sup:`th` configuration in the (2,2) active space. Thus, the
user can easily compare the energy of singlet and triplet microstates.::

  --------------------------------------------------
   Final SA-REKS(2,2) energy:      -4.75910919

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.77753765   1.000000  1.000000  0.00
    OSS   -4.74068073   1.000000  1.000000  0.00
   Trip   -4.77753837   1.000000  1.000000  1.00
  --------------------------------------------------

In this example, for a distorted structure of ethylene, the energy of triplet
microstate is close to that of the PPS state, since the frontier orbitals are
the localised :math:`\pi` orbitals in this system. If one want to know the
gradient of the triplet microstate as well as its relaxed density, these can be
obtained by using ``TargetMicrostate`` keyword in ``REKS`` block. If the value
for this keyword is to ``5``, then the properties will be calculated according
to the index of the microstate. In a (2,2) active space, the 5\ :sup:`th`
microstate indicates a triplet configuration, thus the output shows the
quantities for this microstate. The following results are obtained from the
distorted structure of ethylene molecule::

  --------------------------------------------------
   Gradient Information
  --------------------------------------------------
   5 microstate
      -0.00926010    -0.00000000     0.00000000
       0.00926011     0.00000000    -0.00000000
       0.00096561     0.00066166     0.00000008
       0.00096561    -0.00066166    -0.00000008
      -0.00096562    -0.00000009     0.00066176
      -0.00096562     0.00000009    -0.00066176
  --------------------------------------------------

The gradient is now calculated for the 5\ :sup:`th` microstate. In addition, the
energy of spin contribution in the *detailed.out* file is -0.023, Hartree which
corresponds to the spin constant :math:`W_{pp}` for a carbon atom. In this case
the frontier orbitals consisted of only `p` orbitals of a carbon atom, thus the
energy of the spin contribution mostly consists of interactions between these
orbitals. With this option, one can run the molecular dynamics simulation for
the triplet microstate.

This example shows an input file for calculation of a triplet microstate::

  Reks = SSR22 {
    Energy = {
      Functional = { "PPS" "OSS" }
    }
    TargetState = 1
    TargetMicrostate = 5
    FonMaxIter = 100
    shift = 20.0
    Gradient = ConjugateGradient {
      CGmaxIter = 100
      Tolerance = 1.0E-8
      Preconditioner = Yes
      SaveMemory = Yes
    }
    VerbosityLevel = 1
  }

Note that ``TargetMicrostate`` keyword can be used only with the SA-REKS input
settings discussed above.
