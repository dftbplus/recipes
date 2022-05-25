.. highlight:: none
.. _md-restart:

*****************************
Restarting molecular dynamics
*****************************

To restart a molecular dynamics (MD) calculation, such that the
trajectory is continuous with the earlier calculation:

* Set up a new calculation with the same model (hamiltonian, boundary
  conditions, ...) parameters as the previous one.

* Use the last geometry of the previous run as an initial geometry.

* Copy the velocities (the last 3 columns) of the final geometry in
  the xyz-file of the last run and add them as the initial velocities
  into the dftb_in.hsd file of the new run. (Make sure you set the
  unit of the velocities to be [AA/ps]).

* If you have a temperature profile (see :ref:`md-sim-anneal`), you
  will have to adjust it to start from the right temperature in the
  new run.

Thermostats
-----------

The internal state of the thermostat (if used) can affect the
trajectory of the restart::

  * Berendsen -- no internal state, so the restarted trajectory is
    continuous if restarted, if compared to a longer run from the
    initial conditions.

  * Andersen -- stochastic thermostat, its internal state cannot
    currently be restarted, so the trajectory can be different
    depending on the `RandomSeed` used for the calculation. The
    *average* or ensemble behaviour is consistent however.

  * NoseHoover -- internal state of the thermostat can be specified,
    see the `x`, `v` and `g` values printed in `md.out`. The values
    from the last step should be set in the new input::

      Thermostat = NoseHoover {
        x = { 0.9996930771E+00    0.9964985276E+00    0.9965061982E+00}
        v = {-0.5716463853E-05   -0.3655208165E-04   -0.3640230213E-04}
        g = {-0.7062330098E-07   -0.1886240308E-06   -0.1879742133E-06}
      }


Automating restarts
-------------------

There is a set of `scripts
<https://github.com/korintje/dftbplus_restarter>`_ provided by
Dr. Takuro Hosomi (Department of Applied Chemistry, Graduate School of
Engineering, The University of Tokyo) which can automate this process.

The output can be processed to produce suitable restart file(s), by
running the script in the directory containing the files
`dftb_in.hsd`, `geo_end.xyz` and `charges.bin` (if SCC) from the
initial part of the MD trajectory ::

  python3 ./restart_filemaker.py

This produces a directory called `restart` with the processed new
input file. The script has options that customize its behaviour
(see `restart_filemaker.py -h`). The

A second script, `restart_collector.py`, can be used to collect and
join together the output of several restarted calculations ::

  python3 ./restart_collector.py

This script should be run in the top directory of the calculation, and will
recursively descend to find data ::

  ./geo_end.xyz
    ./restart/geo_end.xyz
      ./restart/restart/geo_end.xyz

and place it into `./collect/` as a combined output file. Again, see
`restart_collector.py -h` for extra options.

Example of restarting
---------------------


.. only :: builder_html or readthedocs

   [Input: `recipes/moleculardynamics/thermalise/`]

Stop a running MD calculation by either

  * pressing `ctrl + c` or in the DFTB+ working directory create a
    stop file ::

      touch stop_driver

    which cleanly halts DFTB+ (remove the file afterwards)

  * waiting for it to terminate after sufficent iterations

Then, in the working directory run the `restart_collector.py`
script. Without any arguments, it will create a directory called
`restart` containing the input to continue the trajectory. Compare the
original `dftb_in.hsd` file and the `restart/dftb_in.hsd`.

You'll notice that the generated `dftb_in.hsd` looks quite different,
but contains the last structure, it's final velocities, any thermostat
state and `ReadInitialCharges = Yes` (as this is an SCC calculation).

This is now ready to restart in the directory (and can be nested to
futher restarts inside this directory).

The resulting sections of trajectories can the be re-assembled using
the `restart_collector.py` script. Run this in the top directory
containing the first calculation in the set, and the resulting
composite data is then placed in the subdirectory `collect/`.
