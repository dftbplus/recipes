.. _sec-linresp:

************
Introduction
************

In linear response TD-DFTB, the transition energies and oscillator
strengths are obtained by solving the so-called Casida equation:

.. math:: \boldsymbol{\Omega} \mathbf{F}_I = \omega_{I}^2 \mathbf{F}_I

:math:`\boldsymbol{\Omega}` is called the response matrix and its
elements depend on the occupied and virtual Kohn-Sham orbitals
(:math:`\phi_{i\sigma}` and :math:`\phi_{a\sigma}`, respectively) and
their energy difference, :math:`\omega_{ia\sigma} =
\epsilon_{a\sigma} - \epsilon_{i\sigma}`:

.. math:: \Omega_{ia\sigma,jb\tau} = \delta_{ij} \delta_{ab}
   \delta_{\sigma \tau} \omega_{jb\tau}^2 + 2 \sqrt{\omega_{ia\sigma}
   \omega_{jb\tau}} K_{ia\sigma, jb\tau}

The matrix :math:`\mathbf{K}` is obtained as the first derivative of
the DFTB Hamiltonian with respect to the ground-state electron density
matrix.

The eigenvalues of this problem are the square of the excitation
energies, while from the eigenvectors, :math:`\mathbf{F}_I`, one can
derive excited state properties such as the excited state total spin
and the oscillator strengths of the transitions.  Specifically, the
latter property can be obtained using:

.. math:: f^I = \frac{2}{3} \sum_{k=1}^3 \left| \sum_{ia\sigma}
   d^k_{ia\sigma} \sqrt{ \omega_{ia\sigma}} F_{ia\sigma}^I \right|^2

where :math:`d^k_{ia\sigma}` are the *k*-th component of the elements
of the dipole matrix.

For close-shell systems, it is possible to unitary-transform the
Casida equation into two independent eigenvalue problems for singlet
and triplet transitions. This reduces the dimension of the equations
and allows us to study excited states with different multiplicities,
independently.

It is important to mention that due to the minimal basis set employed
within DFTB, only valence excited states can be computed with this
formalism.
