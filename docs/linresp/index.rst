.. _linear-response:

###################
Linear response
###################

In this chapter we will go through a couple of recipes for the
computation of excited state properties using linear response
time-dependent DFTB (TD-DFTB).

We can employ TD-DFTB to efficiently calculate vertical and adiabatic
excitation energies, as well as the oscillator strengths of individual
transitions. With this information, it is possible to compute the
absorption spectra of a given system for a defined energy
range. Furthermore, other properties like charges and energy gradients
can be computed for specific excited states.

With DFTB+, we can investigate excited state properties of both
closed- and open-shell molecular compounds. As DFTB+ is currently
limited to finite systems only, the study of periodic structures has
to be worked around by using properly optimised cluster models.  As an
instance of the latter case, we will examine below a recipe for the
computation of the optical spectra of a titania-based cluster system.

Moreover, DFTB+ allows to speed up linear response calculations by
reducing the dimension of the eigenvalue problem. This can be done for
certain systems without incurring on a significant loss of accuracy,
as will be shown in our macromolecule recipe.

Next, we will give a very brief introduction to linear response
TD-DFTB. Feel free to skip this section if you are familiar with the
theory. Afterwards, recipes for some diatomic molecules will be
provided, covering different ground state multiplicity cases. We will
then examine more complex systems, where the real power and efficiency
of TD-DFTB will be exploited.

.. toctree::
   :maxdepth: 1

   introduction.rst
   diatomic.rst
   macromolecule.rst
   no-titania.rst
   rangeseparated.rst
   relax.rst
