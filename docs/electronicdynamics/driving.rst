.. highlight:: none
.. _sec-driving:

************************************************
Driving electronic dynamics with external fields
************************************************

[Input: `recipes/electronicdynamics/driving/`]

In this section we will describe how to drive the Q band of chlorophyll a using
both a continuous and pulsed laser. The first step is to find the transition
dipole moment direction of the Q band. This is done using the
``calc_timeprop_maxpoldir`` script either available after `make install` of
DFTB+, of located in the ``tools/misc`` directory under the ``dftbplus`` source
tree. The invocation of the script is as follows::

  calc_timeprop_maxpoldir -d 20.0 -w 636.0

which produces the output::

  PolarizationDirection = -0.08808129 0.99564018 -0.03069709

these are the three Cartesian components of the transition dipole moment
vector. This vector points in the direction in which a driving laser, if it is
in resonance with the excitation, will have maximum absorption. This vector is
the eigenvector of the polarizability tensor at that energy that has the largest
eigenvalue. It is equivalent to a principal axis of inertia in rigid body
rotation.

Using the ``ElectronDynamics`` block that follows::

  ElectronDynamics = {
    Steps = 60000
    TimeStep [au] = 0.2
    Perturbation = Laser {
      PolarizationDirection = -0.08808129 0.99564018 -0.03069709
      LaserEnergy [eV] = 1.94944
    }
    FieldStrength [V/A] = 0.0001
  }

We resonantly excite the Q band along the direction of its maximum
polarizability. The obtained magnitude of the dipole moment as a function of
time is shown in the following figure:

  .. figure:: ../_figures/elecdynamics/muvst.png
     :height: 40ex
     :align: center
     :alt: Dipole moment as a function of time.

An initial transient, the dipole moment maxima and minima grow in absolute value
as a linear function of time after, confirming that the applied field is in
resonance with the excitation and within the linear response regime. The slope
of this growth is related to the transition dipole of the excitation.

Off-resonant excitation
=======================

[Input: `recipes/electronicdynamics/driving-oot/`]

An *off resonance* excitation at 1.9 eV produces the following result for the
dipole moment:

  .. figure:: ../_figures/elecdynamics/muvst-oot.png
     :height: 40ex
     :align: center
     :alt: Dipole moment as a function of time.
	   
showing characteristic *beats*, the frequency of which are related to the amount
of *detuning*. The amplitude of the dipole moment change caused by the
illumination is also much smaller that when the laser is *in tune* with the
excitation.
