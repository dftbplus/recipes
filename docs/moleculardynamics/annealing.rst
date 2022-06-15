.. highlight:: none
.. _md-sim-anneal:

*******************
Simulated annealing
*******************

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/annealSW/`]

By heating up a system, and then allowing it to slowly cool, stable minima can
be found. If the cooling is sufficiently (logarithmically) slow, then the global
minima will be found, however this is not usually practical. More realistically,
several heating and cooling experiments can find some stable structures, or
investigate the stability of known geometries.

DFTB+ can set a sequence of temperatures during an MD calculation
inside the thermostat input block ::

  Temperature [Kelvin] = TemperatureProfile {
    constant      1    100.0
    linear      499   5500.0
    constant    200   5500.0
    linear      500    100.0
  }

Here, a starting temperature of 1 K is heated over 500 steps up to
5500 K, held at that temperature for 200 steps and then cooled down
again to 100 K in 500 steps.

The first example anneals away a Stone-Wales defect in a graphene sheet. This is
acheived by heating the system up to a high (5000K) temperature where C-C bonds
start to break, holding it at these temperatures for a while and then cooling to
a low temperature. The high temperature is above the usual disociation
temperature of graphene, but since only a relatively short (2.4 ps) time is
computed, using a higher temperature accelerates the annealing.


Annealing to multiple minima
----------------------------

.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/annealV2/`]
   
The annealing process can also be used to find alternative local
minina. Here two vacancies in a graphene sheet are heated. Depending
on the initial conditions, this system anneals to different final
structures. The starting velocities are chosen at random, so depending
on the seed value for the random number generator (several different
cases are given in the input file) a different final defect geometry
is obtained by the same cycle of heating and cooling.
