.. highlight:: none

*******************
Electronic dynamics
*******************

Electron nuclear dynamics within the Ehrenfest *ansatz* is based on the integration of the following equation of motion (EOM) for the reduced one body density matrix: 

:math:`\dot{\rho} = \underbrace{-\mathrm{i} (S^{-1} H \rho - \rho H S^{-1})}_\text{electronic terms} - \underbrace{(S^{-1} D \rho + \rho D^\dagger S^{-1})}_\text{non adiabatic terms}`

As noted above the EOM contains a puerely electronic driving term in which the dynamics may be generated either from an external time dependent driving field in :math:`H` or from the fact that the density matrix might be in a non equilibrium state and therefore does nos conmute with :math:`H`. The second term is part of the Ehrenfest *ansatz* and includes the non adiabatic coupling matrix :math:`D`, the elements of which are defined as :math:`D_{\kappa \nu} = \langle \phi_\kappa | \dot{\phi_\nu} \rangle`, where :math:`| \phi_\kappa \rangle` and :math:`| \phi_\nu \rangle` are a pair of atomic orbitals and, in a localised basis set such as the one used in DFTB+, the time derivative of an orbital with respect to time is a function of the velocity of the nuclei over which the orbital lies. The non adiabatic term introduces a mechanism by which nuclear motion can induce electronic transitions.

The forces used to propagate the nuclei are calculated from the instantaneous expectation value of the force operator adding the corresponding self consistent and repulsive terms. This force comes from a non equilibrium density matrix and is the form in which an excited density matrix may induce nuclear motion.

The EOM for the density matrix is integrated using a Leapfrog scheme in which the density matrix at time :math:`t_{i+i}` is obtained from its value a t time :math:`t_{i-1}` and its derivative at time :math:`t_i`:

:math:`\rho_{i+1}=\rho_{i-1}+2\Delta t \dot{\rho}_i`

The nuclear motion is integrated using a velocity Verlet algorithm. 

The propagation of the electronic and nuclear dynamics may be used for the calculation of absorption spectra, to study the response to constant or pulsed illumination and even the simulation of pump-probe spectroscopy.

All input related to electronic dynamics is located within the *ElectronDynamics* block as part of the *Analysis* block of the input file. Electron dynamics is run after self consistency is achieved, with the possibility of reading previously converged charges for the ground state. Different perturbations can be used, such as A Dirac Delta type pulse, used to obtain the absoprtion spectrum, or a pulse with a Gaussian envelpe as two possible examples.
