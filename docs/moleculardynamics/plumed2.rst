.. highlight:: none

************************************
PLUMED2 integration and metadynamics
************************************

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/plumed2/`]

PLUMED2 is a molecular dynamics biasing and analysis package that has been 
interfaced with DFTB+ to give metadynamics and free energy sampling 
functionality. To use PLUMED2, set the Plumed tag in the Driver block to 
Yes when VerlocityVerlet is selected as the driver. a `plumed.dat` file
must be included in the run directory ::
 
 Driver = VelocityVerlet{
  TimeStep [fs] = 1.0

  Plumed = Yes

  Thermostat = NoseHoover {
    Temperature [Kelvin] = 400
     CouplingStrength [cm^-1] = 3050
 }
 .
 .
 .
 }

The `plumed.dat` file contains information to be read by the PLUMED2
code that will allow it to either bias or analyse the molecular
dynamics on the fly. ::

 DISTANCE ATOMS=4,9 LABEL=d1
 DISTANCE ATOMS=5,9 LABEL=d2

 METAD ...
 LABEL=met ARG=d1,d2 PACE=100 HEIGHT=3 SIGMA=0.01,0.01
 FILE=HILLS
 BIASFACTOR=4 TEMP=400
 ... METAD

 PRINT ARG=d1,d2 STRIDE=100 FILE=plumed_o.dat

 ENDPLUMED

The first 2 lines of the example `plumed.dat` file tell PLUMED2 to define 
the distance between atoms 4 and 9 and the distance between atoms 5 and 9 
and label them `d1` and `d2` respectively. the `METAD` tag tells plumed to
read in parameters for a metadynamics bias calculation with `d1` and `d2` as
input collective variables. This metadynamics simulation places a Gaussian
hills function every 100 md steps. `HEIGHT`, `SIGMA` and `BIASFACTOR` are 
all metadynamics parameters while `TEMP` should correlate with the thermostat
temperature. The Gaussian bias functions are written to the file specified 
by the `FILE` tag. The `PRINT` function tells PLUMED2 to print `d1` and `d2`
every 100 steps to the file specified by the `FILE` tag. Finally, the `ENDPLUMED`
tag indicates to PLUMED2 that input is complete.

More complete information on the functionality of PLUMED2 can be found in the
user manual at ``www.plumed.org/doc-v2.5/user-doc/html/index.html``.
