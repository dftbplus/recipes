.. highlight:: none

************************************************
Driving electronic dynamics with external fields
************************************************

In this chapter we will decribe how to drive the Q band of Chlorophyl a using a continuous and pulsated laser. The first step is to find the transition dipole monet direction of the Q band. This is done using the ``calc_timeprop_maxpoldir`` script l located in the ``utils`` directory under the ``dftbplus`` source tree. The invocation of the script is as follows::

  calc_timeprop_maxpoldir -d 20.0 -w 636.0

which produces the output::

PolarizationDirection = -0.08808129 0.99564018 -0.03069709

this are the three cartesian components of the transition dipole moment vector. This vector points in the direction in which a driving laser in resonance with the excitation will have maximum absorption. This vector is the eigenvector of the polarizability tensor at that energy that has the largest eigenvalue. It is equivalent to a principal axis of inertia in rigid body rotation. 

Using the ``ElectronDynamics`` block that follows::


  ElectronDynamics = {
    Steps = 60000
    TimeStep [au] = 0.2
    Perturbation = Laser {
      PolarizationDirection = -0.08808129 0.99564018 -0.03069709
      LaserEnergy [eV] = 1.94944
      }
  FieldStrength [v/a] = 0.0001
}

We excite the Q band, in resonance, in the direction of maximum polarizability. The obtained dipole moment modulus as a function of time is shown in the following figure:

  .. figure:: ../_figures/elecdynamics/muvst.png
     :height: 40ex
     :align: center
     :alt: Dipole moment as a function of time.

The dipole moment maxima and minuma grow in absolute value as a linear function of time after an initial transient, confirming that the applied field is in resonance with the excitation and within the lienar response regime. The slope of this growth is related to the transition dipole of the excitation.

An *out of tune* excitation at 1.9 eV produces the following result for the dipole moment:

  .. figure:: ../_figures/elecdynamics/muvst-oot.png
     :height: 40ex
     :align: center
     :alt: Dipole moment as a function of time.

showing characteristic *beats*, the frequency of which is related to the amount of *detuning*. The amplitude of the dipole moment change caused by the illumination is also much smaller that when the laser is *in tune* with the excitation.






