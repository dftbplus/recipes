.. highlight:: none

***************************************************
Excitation energies of diatomic molecules
***************************************************

In this section, we show how to compute the first electronic excitation
energies and oscillator strengths of three diatomic molecules, namely NO, N2 and
O2. These simple examples cover three different spin multiplicities of the ground
state: doublet, singlet and triplet, respectively.

Nitrogen molecule
=================

[Input: `recipes/linresp/diatomic/N2`]

First, we examine the N2 case. This is a close-shell diatomic molecule. We are
interested in the low-lying electronic transitions to excited states with both
singlet and triplet multiplicities. The solution of the Casida equation is
performed after a single-point (static) DFTB ground-state calculation for the optimised
geometry. All this is done in a single DFTB run (one single input file). The
*dftb_in.hsd* input file should look like this::

  Geometry = GenFormat {
      <<< "geo_in.gen"
  }

  Driver = {}

  Hamiltonian = DFTB {

      SCC = Yes
      SCCTolerance = 1.0E-10
      MaxAngularMomentum = {
          N = "p"
      }
      SlaterKosterFiles = Type2FileNames {
          Prefix = "../../../slakos/mio-ext/"
          Separator = "-"
          Suffix = ".skf"
      }
      SpinConstants = {
          N = {-0.026} # HOMO Wpp
      }
  }

  ExcitedState {
      Casida {
          NrOfExcitations = 10
          Symmetry = both
          Diagonaliser = Arpack{}
      }
  }

The input shows a standard spin-unpolarised static DFTB calculation,
except for a new block denoted as *ExcitedState* with an embedded subblock
called *Casida*. Once the SCC loop is converged, this block instructs DFTB+ to
build the response matrix according to Casida, using the just-obtained DFTB
orbitals and energies, and diagonalise it.

Apart from this new block, there is another unusual thing in our input. The
Hamiltonian block contains the spin constants, although our calculation is not
spin-polarised. Those are actually not used for the ground state calculation,
but necessary for the computation of the triplet excitations. That is, if we
would be only interested in singlet transitions, we could just omit the spin
constants in the input file.

Let us now take a closer look at the *Casida* block. The keyword *NrOfExcitations*
specifies how many transitions per spin symmetry, or multiplicity, we want
to compute. In *Symmetry* we specify the multiplicity of the transition (either
singlet, triplet or both). The multiplicity of the transition is the difference
between the multiplicities of the excited state and the ground state. In our
example, we are instructing DFTB+ to calculate the first 10 excitations with singlet
symmetry and the first 10 excitations with triplet symmetry. Since the N2 ground
state has multiplicity 1, this means the first 10 singlet-to-singlet and 10
singlet-to-triplet transitions.

Once the calculation is finished (it takes a fraction of a second), the output
file *EXC.DAT* contains the excitation energies and oscillator strengths,
as well as other valuable information. It should look like this::


    w [eV]       Osc.Str.         Transition         Weight      KS [eV]    Sym.

  ================================================================================

     7.321        0.00000000         4   ->     7        0.648       9.014      T
     8.118        0.00000000         5   ->     7        1.000       8.118      T
     8.118        0.00000000         5   ->     6        1.000       8.118      T
     9.014        0.00000000         4   ->     6        0.826       9.014      T
     9.014        0.00000000         3   ->     6        0.715       9.014      T
     9.014        0.00000000         3   ->     7        0.930       9.014      T
    12.062        0.00000000         2   ->     6        0.997      12.062      T
    12.062        0.00000000         2   ->     7        0.997      12.062      T
    22.118        0.00000000         1   ->     6        0.795      22.118      T
    22.118        0.00000000         1   ->     7        0.795      22.118      T
     8.118        0.00000000         5   ->     7        0.804       8.118      S
     8.118        0.00000000         5   ->     6        0.804       8.118      S
     9.014        0.00000000         3   ->     7        0.911       9.014      S
     9.014        0.00000000         4   ->     6        0.969       9.014      S
     9.014        0.00000000         4   ->     7        0.772       9.014      S
    12.062        0.00000000         2   ->     7        0.978      12.062      S
    12.062        0.00000000         2   ->     6        0.978      12.062      S
    12.750        0.80120458         3   ->     6        0.647       9.014      S
    22.118        0.00000000         1   ->     6        0.884      22.118      S
    22.118        0.00000000         1   ->     7        0.884      22.118      S


The triplet transitions are listed first, followed by the singlet ones. They can
be identified by the letter *T* or *S* in the last column.

The first column *w [eV]* is the excited state energy we are looking for, the second one *Osc.Str.* lists the corresponding oscillator strength. The column *Transition* reports the indices of the dominant molecular orbitals involved in the electronic transition. In our example, the singlet state at 12.75 eV features a transition from the occupied Kohn-Sham orbital 3 (HOMO-2) to the virtual orbital 6 (the LUMO). The next column *Weight* indicates the weight of the corresponding singly excited determinant in the CIS expansion of the excited state. Values close to one indicate that the excited state is well described by a single electronic excitation, while small values speak for a collective excitation. Column *KS [eV]* provides the Kohn-Sham transition energy difference  :math:`\omega_{ia\sigma} = \epsilon_{a\sigma} - \epsilon_{i\sigma}` (see above).

Oxygen molecule
=================

[Input: `recipes/linresp/diatomic/O2`]

For the O2 molecule, we will consider its triplet ground state. This is
specified in the input file through the *Hamiltonian/SpinPolarisation* block::

  SpinPolarisation = Colinear {
      UnpairedElectrons = 2
  }

Our excited state block will in this case looks like this::

  ExcitedState {
      Casida {
          NrOfExcitations = 10
          Diagonaliser = Arpack{}
      }
  }

We are instructing DFTB+ to compute the first 10 excitations. Note that since
our system is not closed-shell, we can no longer separate our eigenvalue problem
in two independent singlet and triplet equations, so we have to build and
diagonalise the entire response matrix in this case. But, how do we know the
spin multiplicities of the computed transitions? We get this information from
the last column of the *EXC.DAT* file::

  w [eV]       Osc.Str.         Transition         Weight      KS [eV]    D<S*S>

  ================================================================================

   6.353        0.00000000      5   ->     6        0.829       6.353     0.000
   6.353        0.00000000      4   ->     6        0.787       6.353     0.000
   6.353        0.00000000      5   ->     7        0.722       6.353     0.000
   6.793        0.00000000      3   ->     6        0.993       6.793     0.000
   6.793        0.00000000      3   ->     7        0.993       6.793     0.000
   8.204        0.23976646      4   ->     7        0.617       6.353     0.007
  14.567        0.00000000      2   ->     7        0.989      14.567     0.000
  14.567        0.00000000      2   ->     6        0.989      14.567     0.000
  22.424        0.00000000      6   ->     8        0.800      22.424     0.000
  22.424        0.00000000      7   ->     8        0.800      22.424     0.000

In the last column are the expectation values of the square of the total spin operator for the
transitions. A value of zero means we have a singlet transition (triplet to
triplet). Note that we may have transitions with some spin contamination
(transitions leading to unphysical states). In our next example, we will explore
this in more detail.

Nitric oxide molecule
=====================

[Input: `recipes/linresp/diatomic/NO`]

Finally, we have the NO molecule, with one unpaired electron (doublet ground state). The symmetry of NO leads to degenerate orbitals, which causes problems with the SCC convergence.
We therefore additionally provide a small electronic temperature to ease the ground state computation::

  Filling = Fermi {
        Temperature [K] = 40
  }

NOTE: If your structure features a band gap, it is typically neither necessary nor advisable to set the electronic temperature different from 0 K. The code also works with fractional occupations, but the response matrix will be much larger then necessary and this will cause long calculations. In addition, the memory need increases significantly.



In this case, the first 10 excitations are::

  w [eV]       Osc.Str.         Transition         Weight      KS [eV]    D<S*S>

  ================================================================================

  7.478        0.00209868      4   ->     7        0.556       8.534     1.418
  7.772        0.00000000      5   ->     6        1.000       7.772    -0.000
  7.772        0.00000000      5   ->     7        1.000       7.772    -0.000
  7.793        0.00000000      5   ->     6        0.779       7.793     0.999
  7.793        0.00000000      5   ->     7        0.779       7.793     0.999
  8.534        0.00000000      4   ->     6        0.707       8.534     0.999
  8.534        0.00000000      3   ->     6        0.728       8.534     0.999
  8.636        0.00000000      4   ->     6        0.984       8.636     0.000
  8.636        0.00000000      3   ->     7        0.657       8.636    -0.000
 11.652        0.49971991      3   ->     6        0.600       8.636    -0.597

Let us pay attention to the last column of the *EXC.DAT* file. Contrary to the
previous case, here we obtain large non-zero :math:`\Delta S^2` values. When
:math:`\Delta S^2 = 0`, we are in the presence of a doublet-to-doublet
transition. Likewise, if :math:`\Delta S^2 = 2`, we would have an excitation to
a quadruplet state. Otherwise, we have some extent of spin contamination in
our obtained transitions. The last column should help us determine which
excitations are to be trusted. We can set an arbitrary spin contamination threshold
to establish which transitions we will consider leading to a physical excited state.
In our *NO-TiO2* recipe, we will compute the absorption spectrum of a system where
transitions with a spin contamination beyond an imposed threshold are excluded.
