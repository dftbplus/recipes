####
REKS
####

This chapter discusses REKS (spin-Restricted Ensemble Kohn-Sham) calculations
for ground and low-lying excited states. REKS can treat states with strong
static correlation (`J. Phys. Chem. A 104, 6628
<https://dx.doi.org/10.1021/jp0002289>`_) and also produces the correct spin
symmetry in cases where unrestricted methods suffer from spin contamination.
REKS also allows the calculation of, for example, states with multi-reference
character as well as the vertical excitation energies between states. REKS can
be classified as :ref:`single_state_REKS`, :ref:`SA_REKS` and :ref:`SI_SA_REKS`.

In this chapter we will

* provide some example calculations of basic energy and gradient evaluation,
* discuss additional functionality such as level shifting and spin tuning,
* demonstrate advanced calculations including non-adiabatic couplings between the
  states and relaxed densities which can be used for QM/MM calculations.

Finally, we will discuss how to diagnose and handle some the problems which can
specifically arise in REKS calculations.

.. toctree::
   :maxdepth: 1

   introduction.rst
   advanced.rst
   error-handling.rst
