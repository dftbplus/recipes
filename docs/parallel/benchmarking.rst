Benchmarking and scalability
============================

DFTB+ has an internal timer for various significant parts of its
calculations. This can be enabled by adding the following option to the input::
  
  options ={
    TimingVerbosity = 2
  }

This will activate the timers for the most significant stages of a
calculation. Higher values increase the verbosity of timing, while a value of
`-1` activates all timers.

The typical output of running the code with timers enabled for a serial
calculation will look something like ::
  
  --------------------------------------------------------------------------------
  DFTB+ running times                          cpu [s]             wall clock [s]
  --------------------------------------------------------------------------------
  Pre-SCC initialisation                 +     3.46  (  5.8%)      3.45  (  5.8%)
  SCC                                    +    47.34  ( 79.6%)     47.34  ( 79.6%)
    Diagonalisation                           44.79  ( 75.3%)     44.79  ( 75.3%)
    Density matrix creation                    2.21  (  3.7%)      2.21  (  3.7%)
  Post-SCC processing                    +     8.68  ( 14.6%)      8.68  ( 14.6%)
    Energy-density matrix creation             0.55  (  0.9%)      0.55  (  0.9%)
    Force calculation                          7.74  ( 13.0%)      7.74  ( 13.0%)
    Stress calculation                         0.94  (  1.6%)      0.94  (  1.6%)
  --------------------------------------------------------------------------------
  Missing                                +     0.00  (  0.0%)      0.00  (  0.0%)
  Total                                  =    59.48  (100.0%)     59.47  (100.0%)
  --------------------------------------------------------------------------------

This shows the computing time (cpu) and real time (wall clock) times of major
parts of the calculation.

An alternative method to obtain total times is by using the built-in shell
command `time` when running the DFTB+ binary ::

  time dftb+ > output

In the above example, at the termination of the code timings will be printed::
  
  real    0m59.518s
  user    0m59.430s
  sys     0m0.092s

More advanced timing is possible by using profiling tools such as `gprof`.
  
Examples
--------

Shared memory parallelism
~~~~~~~~~~~~~~~~~~~~~~~~~

The strong scaling for the OpenMP parallel code for some example inputs is shown
here

  .. figure:: ../_figures/parallel/openMP.png
     :height: 40ex
     :align: center
     :alt: Strong scaling for some SiC supercells.

These results were calculated on an old Xeon workstation (E5607 @ 2.27GHz) and
are a self-consistent calculation of forces for different sized SiC supercells
with 1 k-point and `s` and `p` shells of orbitals on the atoms. Timings and
speed-ups come from data produced from the code's reported total wall clock
time.

There are several points to note

  1. The parallel scalability improves for the larger problems, going from ~68%
     to ~79% on moving from 512 atoms to 1728. This is a common feature, that
     larger problems give sufficient material for parallelism to work better.

  2. The gain in throughput for these particular problems is around a factor of
     2 when using 4 processors, and raises to around 3 on 8 processors for the
     largest problem.

  3. From Amdhal's law we can estimate the saturating limits for large numbers
     of processors as ~3.1 and ~4.8 for the smallest and largest problems
     respectively. This implies that there is not much value in using more than
     ~4 processors for the smallest calculation, since this has already gained
     around 2/3 of it theoretical maximum speed up. Adding a 5th or 6th
     processor will only improve performance by ~5% each, so is probably a waste
     of resources. Similarly For the largest calculation in this example, 6
     processors gives around a factor of 3 speed up compared to serial
     operation, but adding 7th processors will only speed the calculation up by
     ~6%.

  4. The experimental data does not lie exactly on the Amdahl curves, this could
     be due to competing processes taking resources (this is a shared memory
     machine with other things running) or the problem may run anomalously well
     for a particular number of processes (in this example 4 processors
     consistently ran slightly better, perhaps due to the cache sizes on this
     machine).

  5. The weak scaling (increasing the number of processors proportional to the
     number of atoms) shows an approximately :math:`O(N^2)` growth in time,
     where serial operation would be :math:`O(N^3)`.

  6. These timings are for this specific hardware and these particular problems,
     so you should test the case you are interested in before deciding on
     parallel resources.

     
Weak scaling from the same data set is shown here
     
  .. figure:: ../_figures/parallel/weakOpenMP.png
     :height: 40ex
     :align: center
     :alt: Weak scaling for some SiC supercells.
	   
