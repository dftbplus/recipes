###################
Electronic dynamics
###################

This chapter discusses DFTB+ calculations where the electrons are allowed to move.
This allows the calculation of electronic absobption spectra when the perturbing field
is a Dirac delta pulse. More general external time dependent fields are allowed, 
including arbitrary eliptical polarization, and pulses with various envelope functions.

Ehrenfest dynamics can also be performed. In this methods nuclei move driven
by the instantaneous expectation value of the force at each time for the moving
electron density.

In this chapter we will provide example calculations of electronic spectra, driving
of the electronic dynamics using external pulsed and continuous fields and 
driving of electronic motion via electronic excitation.

.. toctree::
   :maxdepth: 1

   introduction.rst
   spectra.rst
   driving.rst
   ehrenfest.rst
