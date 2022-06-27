.. highlight:: none

***********************************************************************
Tutorial first steps: Absorption spectra, analysis and charge transfer
***********************************************************************


This is the tutorial presented in the DFTB+ school 2022 in Daresbury. The idea is to get used
to with the real-time TDDFTB method learning the basics of absorption spectra calculations and 
photoinduced processes like charge transfer under irradiation.

Spectra and analysis
====================

We will first calculate the absorption spectra of two different molecules and
analyse them using the tools provided by DFTB+.

Calculation of the absorption spectrum of carbazole
---------------------------------------------------

We will calculate the absorption spectrum of the carbazole molecule

1. Take a look at the input coordinates *coords.gen*. The *.gen* format
   is the one used for DFTB+ code. In order to visualize the molecule,
   you can use the ``gen2xyz`` script, provided in the installation of the 
   DFTB+ code, doing::
     
     gen2xyz coords.gen

   This will generate a *coords.xyz* that you can open with VMD, Avogadro or
   any other molecular visualization software of your choice.

   .. figure:: ../_figures/elecdynamics/tutorial/carbazole.png
      :height: 30ex
      :align: center
      :alt: Carbazole molecule.

2. Open the *dftb_in.hsd_spec* file. This is a template for the calculation
   of the absorption spectrum.

   - Take a look at the ``ElectronDynamics`` block at the end of the file:: 
    
      ElectronDynamics = {
          Steps =                         #define the time window
          TimeStep [au] =                 #resolution of the spectrum
          Perturbation = Kick {           #must be a kick (Dirac delta)
              PolarizationDirection =     #desired direction/s
          }
          FieldStrength [v/a] =           #Field strength of the perturbation
      }

   - The input variables to be considered for the calculation of the spectrum are four:

     * ``Steps`` (integer): the number of steps of the dynamics. It defines the time window and,
       consequently, the energy window of the spectrum. The longer the dynamic, the lower the
       energies that can be reached in the spectrum are (also afected by the timestep, of course).
       Here, we will use 10000 steps.
     * ``TimeStep`` (float): the time step in time units. Usually 0.2 a.u (0.0048 fs). 
       It defines the resolution of the spectrum. The smaller the time step, 
       the higher the resolution within the time window.
     * ``Perturbation``: In this case we need a kick perturbation (Dirac delta) and we need to 
       specify the ``PolarizationDirection`` that could be *X*, *Y*, *Z* (if we are interested in 
       one particular direction) or *all* if we want to calculate the whole absorption spectrum.
       Set it as *all*.
     * ``FieldStrength`` (float): the field strength of the perturbation applied. For the
       calculation of the absorption spectrum, it must be within the linear response regime,
       i.e. usually 0.001 :math:`V/\AA`.

   - Complete the template file, copy it to *dftb_in.hsd* and run the calculation.

3. Once the dynamics ended, we will have 3 components of the dynamical dipole moment 
   (*mux.dat*, *muy.dat*, *muz.dat*). We need to Fourier-transform these dipole components
   in order to obtain the absorption spectrum of the molecule. To do this, we will use the
   ``calc_timeprop_spectrum`` tool available after installation of DFTB+ under: 
   *dftbplus/tools/misc/*. In the folder
   where you have the dipole files just type::

    calc_timeprop_spectrum -d 4 -f 0.001

   where the option -d is for the damping constant (in fs) applied to the dipole moment before transformation.
   The option -f stands for the field strength (in V/AA) of the perturbation applied during dynamic.

4. After running the script you will find two new files: *spec-nm.dat* and *spec-ev.dat* which are
   the absorption spectra in nm and eV, respectively. Plot the spectrum file with the plotting tool
   of your desire and look at the lower energy transitions. You should then see
   an absorption spectrum similar to:

   .. figure:: ../_figures/elecdynamics/tutorial/spec-nm-carbazole.png
      :width: 60%
      :align: center
      :alt: Absorption spectrum of carbazole molecule

      Absorption spectrum of carbazole molecule

5. Change the damping constant for a higher value, recalculate the specctrum and plot both spectra
   together. Which is the effect of the damping time in the spectrum?
   Here it is an example of the same spectrum obtained before, calculated with
   different values of the damping constant.

   .. figure:: ../_figures/elecdynamics/tutorial/specs-comparison-damp.png
      :width: 60%
      :align: center
      :alt: Influence of the damping constant value ``d`` in the absorption spectrum.

      Influence of the damping constant value ``d`` in the absorption spectrum.

Analysis of the absorption spectrum of carbazole
------------------------------------------------

We will consider a laser perturbation in tune with the lowest energy
transition of the molecule in order to study the photodynamic
process of absorption in this transition. In order to do this, we
need to know the energy of the lowest energy transition of
the molecule (look for it in the spectrum plotted in the previous calculation)
and calculate the direction of maximal polarization of the transition.

1. Open the *dftb_in.hsd_laser* file. This is a template for the calculation
   of a laser perturbation.

   - Take a look at the ``ElectronDynamics`` block at the end of the file:: 
     
      ElectronDynamics = {
         Steps = 10000
         TimeStep [au] = 0.2
         Perturbation = Laser{              # Laser type perturbation
            LaserEnergy [nm] =              # energy of interest
            PolarizationDirection =         # calculate with calc_timeprop_maxpoldir
         }
         FieldStrength [v/a] = 0.001
         Populations = Yes                  # to write populations during dynamic
      }

     Now, the ``Perturbation`` type is a ``laser`` (and not anymore a ``kick``)
     and we need to specify two parameters:
         
         * ``LaserEnergy`` (float): the energy of the applied laser that may be
           the transition energy of interest. This value must be in energy units
           like eV but also nm is possible.
         * ``PolarizationDirection`` (vector): in the case of a laser, the 
           ``PolarizationDirection`` is 3-cartesian components vector in which the 
           laser will be applied. 

     Also note that we turned on the ``Populations`` flag in order to write
     the occupations during the dynamics.   

2. To complete the input template for the laser, we need to provide
   the ``LaserEnergy`` and the ``PolarizationDirection`` of the laser. Based on 
   our previous calculated spectrum, calculate the direction of maximal 
   polarization of the lowest energy transition of the molecule.

   - Help: use the tool ``calc_timeprop_maxpoldir`` already available in
     your installation (under: *dftbplus/tools/misc/*). To know how this
     tool work the user can just type::

      calc_timeprop_maxpoldir -h

   - Along which axes is the direction vector? How is this explained?
  
     - Hint: try to visualize the molecule and see how it is oriented with respect
       to the cartesian axes.
   
   + Solution: If you choose the lower energy transition of carbazole you may do::
      
      calc_timeprop_maxpoldir -10 -w 326

     and you will obtain the following transition dipole vector::
      
      PolarizationDirection = 0.99999221 0.00101174 -0.00381496

     which is essentially paralel to the *X* cartesian direction


3. Prepare the input for the dynamics under a continuous laser perturbation.
   Use the energy transition obtained from the spectrum as the ``LaserEnergy``
   and the vector obtained above as the ``PolarizationDirection`` of the 
   laser.

   - Why we should use this direction instead of any other?

4. After the dynamics, take a look at the *mu.dat* file.

   - Is the dipole moment increasing linearly?

5. Take a look at the *molpopul.dat*
   generated. This file contains the populations projected on the GS orbitals during the dynamics.

   - Which orbitals are involved in the transition?
     Help: you can plot the *molpopul.dat* file using xmgrace::

      xmgrace -nxy molpopul.dat

     Look at the populations at y=2 (occupied orbitals in the GS basis) and find
     which curves are decreasing during the dynamic. These are the orbitals
     being depopulated.
     Look at the populations at y=0 (unoccopied orbitals in the GS basis) and find
     which curves are increasing during the dynamics. These are the orbitals
     being populated.

6. Let's generate those orbitals using ``waveplot``

  - Look at the *waveplot_in.hsd_* template input file for waveplot:

    - Which files are needed?
    - In which orbitals are we interested?

  - After editing and completing this file, just rename it to *waveplot_in.hsd* and run
    ``waveplot`` using your current installed version that probably is at::
       
       $HOME/dftbplus/_build/app/waveplot/waveplot

  - After running waveplot, a number of files would be generated starting with "wp-1-1".

7. Let's plot these orbitals:

   - Open the cube files that correspond to the HOMO and LUMO and plot them as an isosurface.
     (For a tutorial on the [Basics of VMD](https://www.ks.uiuc.edu/Training/SumSchool/materials/sources/tutorials/01-vmd-tutorial/html/node2.html) and/or plotting an [isosurface](https://www.ks.uiuc.edu/Research/vmd/current/ug/node77.html) method please refer to the links.)

Here we show a figure with the Populations obtained from the laser dynamics
and the orbitals involved in the transition. You should get something 
similar in your calculations:

.. figure:: ../_figures/elecdynamics/tutorial/molpopul-carbazole.png
   :width: 60%
   :align: center
   :alt: molpopul1-carbazole

   (left)Populations vs time for the laser dynamics. (right) Orbitals involved
   in the lower energy transition of the carbazole molecule.
