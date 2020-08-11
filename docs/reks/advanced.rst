
============================================
Advanced calculations
============================================

************************
Nonadiabatic couplings
************************

The nonadiabatic coupling vectors between the states can be calculated with SSR method.
By setting ``NonAdiabaticCoupling = Yes`` in ``REKS`` block, the user can easily
calculate the nonadiabatic coupling vectors. When this option is turned on, SSR shows
the gradients for all states and possible coupling vectors::

  --------------------------------------------------
   Gradient Information
  --------------------------------------------------
   1 st state (SSR)
       0.00008470    -0.00007273    -0.00001868
      -0.00008072     0.00008276     0.00002092
      -0.00001026     0.00004220    -0.00002731
       0.00004647    -0.00001374     0.00003727
       0.00001848     0.00000383     0.00002557
      -0.00005867    -0.00004232    -0.00003777

   2 st state (SSR)
       0.28581286    -0.00360737     0.00128038
      -0.28580887     0.00361741    -0.00128220
      -0.00001066     0.00004477     0.00015391
       0.00004763    -0.00001656    -0.00016019
       0.00001888     0.00000173    -0.00015362
      -0.00005984    -0.00003999     0.00016172
  --------------------------------------------------
   Coupling Information
  --------------------------------------------------
   between  1 and  2 states
  --------------------------------------------------
   g vector - difference gradient
      -0.14286408     0.00176732    -0.00064953
       0.14286407    -0.00176733     0.00065156
       0.00000020    -0.00000129    -0.00009061
      -0.00000058     0.00000141     0.00009873
      -0.00000020     0.00000105     0.00008960
       0.00000059    -0.00000117    -0.00009975

   h vector - derivative coupling
      -0.00052064     0.00000662    -0.00000429
      -0.00052061     0.00000625    -0.00000486
       0.00025848    -0.00015640     0.00000578
       0.00026214     0.00014996    -0.00000135
       0.00026234     0.00015000     0.00000197
       0.00025828    -0.00015643     0.00000276

   G vector - GDV
      -0.14286408     0.00176732    -0.00064953
       0.14286407    -0.00176733     0.00065156
       0.00000020    -0.00000129    -0.00009061
      -0.00000058     0.00000141     0.00009873
      -0.00000020     0.00000105     0.00008960
       0.00000059    -0.00000117    -0.00009975

   H vector - DCV - non-adiabatic coupling
      -0.00205212     0.00002609    -0.00001692
      -0.00205203     0.00002462    -0.00001917
       0.00101884    -0.00061646     0.00002279
       0.00103324     0.00059107    -0.00000534
       0.00103403     0.00059125     0.00000776
       0.00101804    -0.00061657     0.00001087
  --------------------------------------------------

These gradients and coupling vectors are obtained with planar structure of a ethylene molecule.
The g vector is defined as gradient difference vector, thus it can be calculated from the
difference of SA-REKS gradients. Similarly to this, G vector is calculated from the difference
of SSR gradients. The h vector is defined as coupling gradient, so it can be simply calculated
from the gradient of state-interaction terms. The H vector is defined as derivative coupling
vectors, thus its norm increases as the energy gap becomes smaller.

The g and h vectors can be regarded as the vectors defined through diabatic states, and G and
H vectors are defined through the adiabatic (SSR) states. In general, the nonadiabatic coupling
vectors can be used for Ehrenfest or surface hopping molecular dynamics, while the g and h vectors
can be used for minimum energy conical intersction (MECI) optimisation.

************************
Relaxed Density
************************

The relaxed density is calculated according to ``TargetState`` when the user sets ``RelaxedDensity``
to ``Yes``. The calculation of relaxed density requires the information about gradient, thus it
can be calculated when the input includes calculation of gradient. When this option is turned on,
the relaxed FONs are given in the bottom of standard output and *relaxed_charge.dat* file is
generated. It includes total charge as well as the mulliken charges of each atom for target state.
This example shows relaxed charges for ground state::

  total charge:     -0.00000000 (e)

  relaxed S0 atomic charge (e)
     atom        charge
        1      -0.17586135
        2      -0.17586344
        3       0.08793393
        4       0.08792776
        5       0.08797812
        6       0.08788497

If one want to exploit SSR with QM/MM approach, ``RelaxedDensity`` option should be used to
obtain the gradient of external point charges. Then, the gradient of external point charges
for target state are given in *detailed.out* file. Therefore, this option can be used to run
Ehrenfest or surface hopping dynamics with QM/MM approach.

************************
Spin tuning constants
************************

DFTB/SSR method can well describe equilibrium geometries as well as vertical excitation energies
compared with SSR/wPBEh result. (See the paper `JCTC, 2019, 15, 3021-3032.
<https://pubs.acs.org/doi/10.1021/acs.jctc.9b00132>`_) However, the behaviours at MECI points
does not sometimes match those obtained with SSR/wPBEh. For example, n/:math:`\pi^*` type MECI
geometry of ethylene or methyliminium molecule cannot be located with DFTB/SSR, and different
calculation results mainly originate from an incorrect description of the relative stability of
the PPS and OSS states. Their relative stability depends on the splitting between the open-shell
singlet microstates and the triplet microstates in the PPS and OSS energies.

In principle, DFTB/SSR employs spin-polarized DFTB formalism, in which the spin-polarization
contribution is obtained from the second-order expansion of the magnetization density with 
respect to zero magnetization electronic structure. At the n/:math:`\pi^*` type MECI geometries,
both frontier orbitals are located on the same atom. In such a case, the second-order expansion
of magnetization may not be suitable for the triplet microstates, as the spin density becomes
too large. As a simple solution, the stability between the PPS and OSS states can be adjusted
by scaling the atomic spin constants. For most molecules the FONs for PPS state become :math:`n_a`
~ 2.0 and :math:`n_b` ~ 0.0, hence the energy of the PPS state is determined by 1st microstate
alone and it is the energy of the OSS state that depends the atomic spin constants. If the user
run the test calculation included in $DFTB/test/prog/dftb+/reks/PSB3_2SSR_rangesep_tuning directory,
then followng results can be given in the standard output::

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

It shows the energies at MECI point of PSB3 molecule, thus the nonadiabatic coupling vectors show
large elements near center C=C bond. Similary to this, the atomic spin constants can be modified
by using ``SpinTuning`` keyword in ``REKS`` block as follows::

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

************************
Microstate calculation
************************

Obviously SSR method treats only singlet state such as PPS or OSS state. If one want to compare
the energy of singlet and triplet states, SSR provides the energy of triplet configuration as
an alternative which corresponds to 5th or 6th configuration in (2,2) active space. Thus, the user
can easily compare the energy of singlet states and triplet microstate.::

  --------------------------------------------------
   Final SA-REKS(2,2) energy:      -4.74617355

   State     Energy      FON(1)    FON(2)   Spin
    PPS   -4.76360999   1.244973  0.755027  0.00
    OSS   -4.72873711   1.000000  1.000000  0.00
   Trip   -4.76302402   1.000000  1.000000  1.00
  --------------------------------------------------

In this example showing the resulting energy of distorted structure of ethylene, the energy of
triplet microstate is almost similar with that of PPS state since the frontier orbitals are the
localized :math:`\pi` orbitals in this system. If one want to know the gradient of triplet
microstate as well as relaxed density, then they can be obtained by using ``TargetMicrostate``
keyword in ``REKS`` block. If the value for this keyword sets to ``5``, then the properties will
be calculated according to the index of microstate. In (2,2) active space, 5th microstate indicates
triplet configuration, thus the output shows the quantities for this microstate. The following
results are obtained from the distorted structure of ethylene molecule::

  --------------------------------------------------
   Gradient Information
  --------------------------------------------------
   5 microstate
       0.28582859    -0.00360772     0.00127630
      -0.28582459     0.00361760    -0.00127747
      -0.00001027     0.00004031    -0.00020014
       0.00004510    -0.00001196     0.00019623
       0.00001850     0.00000643     0.00020011
      -0.00005733    -0.00004467    -0.00019502
  --------------------------------------------------

The gradient is now calculated for 5th microstate. In addition, the energy of spin contribution
in *detailed.out* file is -0.023 Hartree which corresponds to the spin constant :math:`W_{pp}` for
carbon atom. In this case the frontier orbitals are consisted of only p orbitals of carbon atom,
thus the energy of spin contribution is mostly consisted of interactions between p orbitals of
carbon atom. With this option, one can run the molecular dynamics simulation for triplet microstate.

This example shows input file for calculation of triplet microstate::

  Reks = SSR22 {
    Energy = {
      Functional = { "PPS" "OSS" }
    }
    TargetState = 1
    TargetMicrostate = 5
    FonMaxIter = 30
    shift = 0.3
    Gradient = ConjugateGradient {
      CGmaxIter = 100
      Tolerance = 1.0E-8
      Preconditioner = Yes
      SaveMemory = Yes
    }
    RelaxedDensity = Yes
    VerbosityLevel = 1
  }

Note that ``TargetMicrostate`` keyword can be used with only SA-REKS input settings as above.


