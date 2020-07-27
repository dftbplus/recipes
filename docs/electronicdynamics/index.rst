###################
Electronic dynamics
###################

This chapter discusses DFTB+ calculations where electrons are allowed to move
(evolve in time).  This allows the calculation, for example, of electronic
absorption spectra if the electrons are initially perturbed from the
ground-state by a Dirac delta pulse. More general external time dependent fields
can also be applied, including arbitrary elliptical polarisation and pulses with
various envelope functions.

Ehrenfest dynamics can also be performed, where both the electrons and ions
move. In this methods nuclei are driven by the instantaneous expectation value
of the force at each time for the moving electron density.

In this chapter we will provide example calculations of electronic spectra,
driving of the electronic dynamics using external pulsed and continuous fields
and driving of nuclear motion via electronic excitation.

.. toctree::
   :maxdepth: 1

   introduction.rst
   spectra.rst
   driving.rst
   ehrenfest.rst
