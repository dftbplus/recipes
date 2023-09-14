**********************
Extended tight binding
**********************

[Input: `recipes/xtb/`]

This chapter will investigate some basic usage of the extended tight binding (xTB) methods.\ :cite:`bannwarth2020`


Calculation setup
=================

To make use of the extended tight binding Hamiltonian, the ``xTB`` group is set for the ``Hamiltonian`` group in the ``dftb_in.hsd``.

.. tab-set::

   .. tab-item:: GFN1-xTB\ :cite:`grimme2017-1`

      .. literalinclude:: ../_archives/recipes/xtb/gfn1/dftb_in.hsd
         :caption: dftb_in.hsd
         :language: shell
         :emphasize-lines: 6

   .. tab-item:: GFN2-xTB\ :cite:`bannwarth2019`

      .. literalinclude:: ../_archives/recipes/xtb/gfn2/dftb_in.hsd
         :caption: dftb_in.hsd
         :language: shell
         :emphasize-lines: 6

   .. tab-item:: IPEA1-xTB\ :cite:`asgeirsson2017`

      .. literalinclude:: ../_archives/recipes/xtb/ipea1/dftb_in.hsd
         :caption: dftb_in.hsd
         :language: shell
         :emphasize-lines: 6

An example system with several different elements can be this molecule of the mindless benchmark set from the GMTKN55:

.. literalinclude:: ../_archives/recipes/xtb/structs/mindless42.xyz
   :caption: struc.xyz

For GFN2-xTB we will find the following printout for the initialization.

.. code-block:: text
   :emphasize-lines: 21-22,30,32-39

   ...

   Starting initialization...
   --------------------------------------------------------------------------------
   OpenMP threads:              4
   Chosen random seed:          1344338798
   Mode:                        Static calculation
   Self consistent charges:     Yes
   SCC-tolerance:                 0.100000E-04
   Max. scc iterations:                    100
   Shell resolved Hubbard:      Yes
   Spin polarisation:           No
   Nr. of up electrons:            28.500000
   Nr. of down electrons:          28.500000
   Periodic boundaries:         No
   Electronic solver:           Relatively robust
   Mixer:                       Broyden mixer
   Mixing parameter:                  0.200000
   Maximal SCC-cycles:                     100
   Nr. of chrg. vec. in memory:            100
   tblite library version:      0.2.1
   -> parametrization:          GFN2-xTB
   -> repulsion:                Yes
   -> dispersion:               Yes
   -> halogen bonding:          No
   -> electrostatics:           Yes
      -> isotropic:             Yes
      -> anisotropic:           Yes
      -> third-order:           Yes
   Electronic temperature:              0.950045E-03 H      0.258520E-01 eV
   Initial charges:             Set automatically (system chrg:   0.000E+00)
   Included shells:             B:  s, p
                                H:  s
                                O:  s, p
                                N:  s, p
                                F:  s, p
                                P:  s, p, d
                                C:  s, p
                               Al:  s, p, d
   Extra options:
                                Mulliken analysis
                                Force calculation
   Force type                   original

   --------------------------------------------------------------------------------

   ...

Three notable points can be found in the output:

1. the version of the library used for evaluating xTB, as well as details on the parametrization of the xTB Hamiltonian and the included interactions

2. the electronic temperature defaulting to 300K with the xTB Hamiltonian

3. the availability of parameters for all the contained elements


.. admonition:: Exercise
   :class: info

   Perform a few trial calculation with the available xTB Hamiltonians.
   Try to check a few larger systems as well.

   .. dropdown:: Taxol molecule

      .. literalinclude:: ../_archives/recipes/xtb/structs/taxol.xyz

   Do you notice a difference in computational cost of the xTB methods?
   Is this difference expected and how can it be explained?


Using parameter files
=====================

For the silicon element a reparametrization of the GFN1-xTB Hamiltonian was proposed, which is not directly available as ``Method`` keyword.
Instead a ``ParameterFile`` can be provided.
To obtain the parameter file for the GFN1(Si)-xTB method,\ :cite:`komissarov2021` we start by creating the GFN1-xTB parameter file from the internal parametrization storage

.. code-block:: text

   ❯ tblite param --method gfn1 --output gfn1-si-xtb.toml
   GFN1-xTB
   S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput., 2017, 13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118
   [Info] Parameter file written to 'gfn1-si-xtb.toml'

Checking the publication for the silicon reparametrization, we find the following block should replace the original silicon parameters

.. literalinclude:: ../_archives/recipes/xtb/si/gfn1-si-xtb.toml
   :caption: GFN1(Si)-xTB silicon parameters
   :language: toml
   :lines: 1157-1175

Furthermore, a new element-pair scaling is added for silicon–oxygen pairs under ``hamiltonian.xtb.kpair``, replacing the default scaling of 1.0 by

.. code-block:: toml

   # ...
   [hamiltonian.xtb.kpair]
   # ...
   Si-O = 0.9689382911964309
   # ...

Finally, it would be diligent to properly update the ``meta`` entry to correctly identify the parametrization

.. code-block:: toml

   [meta]
   name = "GFN1(Si)-xTB"
   reference = """S. Grimme, C. Bannwarth, P. Shushkov, J. Chem. Theory Comput., 2017, 13, 1989-2009. DOI: 10.1021/acs.jctc.7b00118
   L. Komissarov, Toon Verstraelen, J. Chem. Inf. Mod. ASAP. DOI: 10.1021/acs.jcim.1c01170"""

The modified parameter file can now be used in DFTB+ with

.. literalinclude:: ../_archives/recipes/xtb/si/dftb_in.hsd
   :caption: dftb_in.hsd
   :language: shell
   :emphasize-lines: 6

.. tip::

   The specification of the parameter file format can be found in the `tblite documentation <https://tblite.readthedocs.io>`__.


Summary
=======

.. admonition:: You learned...
   :class: important

   - to setup a calculation with the extended tight binding Hamiltonian
   - perform single point evaluations on various element combinations with xTB
   - export parameter files and manipulate the parametrization data
