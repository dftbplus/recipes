.. _sec-parallel-tuning:

Tuning distributed calculations
===============================

The :ref:`previous section <sec-parallel-benchmarking>` covered choosing *how
many* processors to use. This section covers the settings that determine how
well those processors are actually used. For a large self-consistent calculation on
a distributed memory machine, the difference between a default invocation and a
tuned one can be a factor of two, without changing the answer.

All of the numbers quoted below were measured for a single test case: a
16,000 atom SiC supercell (64,000 basis functions, closed shell, single
k-point, `s` and `p` shells), running on four nodes of an InfiniBand connected
cluster, each node holding two 28 core Intel Xeon Platinum 8280 (Cascade Lake)
processors, with Intel MPI and MKL. They are given to show the *size* of the
effects, not as recommended values -- every one of them depends on the machine,
the libraries and the problem. The method of finding them is the part that
transfers to your system.


Start by finding out where the time goes
----------------------------------------

Before tuning anything, enable the internal timers (see
:ref:`sec-parallel-benchmarking`)::

  Options {
    TimingVerbosity = 2
  }

For large self-consistent problems, the dense eigensolver usually dominates::

  SCC                                    +   285.32  ( 80.6%)
    Diagonalisation                          283.02  ( 79.9%)
    Density matrix creation                   42.70  ( 12.1%)

With a breakdown like this, effort spent on the eigensolver and on the
communication underneath it is worthwhile, and effort spent elsewhere is not.
If instead the time is dominated by force evaluation or by the Hamiltonian
build, the settings below will not help you much.


Choosing a solver
-----------------

The solver is set in the Hamiltonian block::

  Electronic {
    Solver = ELPA {
      Mode = 2
    }
  }

For large dense problems on distributed memory machines, the two stage ELPA
solver (`Mode = 2`) was the fastest of the available dense solvers in our
tests, ahead of the ScaLAPACK based `DivideAndConquer`, `QR` and
`RelativelyRobust` options. Depending on how DFTB+ was compiled, ELPA is
reached either directly or through the ELSI library.

Two cautions:

#. **Linear scaling solvers are not automatically faster.** DFTB+ can also use
   the density matrix purification solver `NTPoly` and the pole expansion
   solver `PEXSI` through ELSI. Both have better asymptotic scaling than
   diagonalisation, but they only win above a crossover system size. For the
   16,000 atom example above they were roughly 11 and 15 times *slower* than
   ELPA respectively. If you want to use them, measure the crossover for your
   own class of system rather than assuming it has been passed.

#. **Autotuning has to be paid for.** ELPA can tune its internal kernel choice
   at runtime, which DFTB+ exposes as::

     Solver = ELPA {
       Mode = 2
       Autotune = Yes
     }

   The tuning takes several diagonalisations to converge, so it only pays off
   if there are many of them. For a self-consistent calculation converging in
   about ten iterations it cost us around 21% in total wall clock time. It is
   worth testing for long molecular dynamics runs, where the same matrix shape
   is diagonalised thousands of times.


Distributing ranks and threads
------------------------------

DFTB+ can be built with MPI and OpenMP together. Using both means deciding how
to split the cores of a node between MPI processes and threads per process.
Hybrid operation has to be requested explicitly::

  Parallel {
    UseOmpThreads = Yes
  }

Without this keyword an MPI parallel binary will stop if `OMP_NUM_THREADS` is
greater than one, which is a deliberate safety net against accidentally
oversubscribing a machine.

For our test case, filling the 56 cores of each node as 28 MPI ranks with 2
OpenMP threads each was consistently faster than 56 single threaded ranks --
by about 8% before the other settings in this section were applied, and about
4% after. The reason is that the dense eigensolver's communication volume grows
with the number of MPI ranks in the process grid, while the shared memory
threads inside a rank cost nothing to communicate with. Pushing further in that
direction is not automatically better either: at some point too few ranks
leaves the BLAS unable to use the cores efficiently. Two to four threads per
rank is a reasonable range to scan.

Whatever split you choose, the total of ranks times threads per node should
equal the number of physical cores. Oversubscribing is much more damaging here
than leaving a core idle.


Process and thread pinning
--------------------------

On a multi-socket node, an MPI rank whose threads drift between sockets will
repeatedly fetch its data across the inter-socket link. Pinning each rank to a
fixed set of cores avoids this. The names of the controlling variables are
specific to the MPI and OpenMP runtimes; for Intel MPI and the Intel OpenMP
runtime::

  export OMP_NUM_THREADS=2
  export I_MPI_PIN_DOMAIN=omp
  export I_MPI_PIN_ORDER=compact
  export KMP_AFFINITY=granularity=core,compact,1

`I_MPI_PIN_DOMAIN=omp` sizes each rank's domain from `OMP_NUM_THREADS`, and
`KMP_AFFINITY` keeps the threads of a rank on neighbouring cores. The
equivalent for OpenMPI and GNU OpenMP is a `--map-by` specification together
with `OMP_PROC_BIND` and `OMP_PLACES`.

Two traps are worth pointing out, because both are silent:

#. **Export the thread count before launching, not inline.** Writing
   ``OMP_NUM_THREADS=2 mpiexec ...`` sets the variable for the child process
   only. Some process managers read it in the launcher when deciding how large
   a pinning domain each rank should get, and if it is not set there yet, every
   rank is pinned to a single core while the program still believes it has two
   threads. The calculation runs, produces correct results and is much slower.
   Use ``export`` on a separate line.

#. **Asynchronous progress threads count against the core budget.** Options
   such as `I_MPI_ASYNC_PROGRESS=1` spawn an extra helper thread *per MPI
   rank*, not per node. If the ranks and their compute threads already fill
   every core, the helper threads have nowhere to run. In our tests, every
   layout that used all cores of a node deadlocked in the first
   self-consistent iteration with this enabled. If you want to try it, leave
   a core per rank free for it.


Block size of the distributed matrices
--------------------------------------

Dense distributed matrices are laid out in a block-cyclic distribution, whose
block size DFTB+ exposes as::

  Parallel {
    Blacs {
      BlockSize = 48
    }
  }

The default is 32. Smaller blocks spread the work more evenly over the process
grid, larger blocks give the local BLAS longer runs of contiguous data. The
optimum depends on the matrix size, the process grid shape and the cache sizes,
so it is worth a short scan over, say, 16, 32, 48, 64 and 128 for a production
sized problem. Moving from 32 to 48 gained about 2% for our test case.

This is a pure input file change, so it is one of the cheapest things to test:
no recompilation is involved.


Kernel selection inside ELPA
----------------------------

ELPA implements its inner kernels several times over for different vector
instruction sets, and picks one at runtime. The choice is made by ELPA itself,
not by DFTB+, and can be overridden through ELPA's own environment variables::

  export ELPA_FORCE_real_kernel=ELPA_2STAGE_REAL_AVX512_BLOCK6

The full enumeration name is required, as shown; short forms are not accepted.
On our AVX-512 capable hardware, forcing this kernel was 3.3% faster overall
than the kernel ELPA selected by itself, purely from the choice of the inner
kernel of the tridiagonalisation back-transformation. The corresponding
variable for complex matrices is `ELPA_FORCE_complex_kernel`.

Which kernel wins is a property of the processor and of the matrix and block
sizes, and the available names depend on how ELPA was compiled. Listing them
for your installation and scanning the plausible ones once, for a given machine
and problem class, is a small investment. If ELPA was built without AVX-512
support, forcing an AVX-512 kernel will simply fail, which is a useful way of
finding out what the library actually contains.


Collective algorithms in the MPI library
----------------------------------------

The eigensolver's communication is dominated by a few collective operations.
MPI libraries implement each collective with several algorithms and choose
between them using built-in heuristics based on message size and process count.
These heuristics are not always right for a given fabric. With Intel MPI the
choice can be forced::

  export I_MPI_ADJUST_ALLREDUCE=8
  export I_MPI_ADJUST_BCAST=9
  export I_MPI_ADJUST_REDUCE=1

Choosing algorithms this way was worth about 2% for our test case. The
numbering is specific to Intel MPI; OpenMPI offers similar control through its
`coll` framework MCA parameters. This is worth a scan if profiling shows a
large fraction of time inside collectives, and is worth nothing otherwise.

Also check that the library is using the fastest available transport rather
than falling back to a generic one, for example with `I_MPI_FABRICS` and
`FI_PROVIDER` for Intel MPI over InfiniBand. A wrong fabric choice is a much
larger effect than any of the tuning above.


A method for scanning settings
------------------------------

The settings above interact, and some of them silently change the calculation
rather than merely speeding it up. A workable procedure is:

#. Establish a baseline and measure its run-to-run scatter by repeating the
   same calculation, unchanged, at least three times. On a busy cluster this
   scatter can easily be a few percent, and any "improvement" smaller than it
   is not an improvement.

#. Change one setting at a time from the baseline, rather than accumulating
   changes. Interactions do exist, but they are much easier to interpret once
   the individual effects are known.

#. For every candidate setting, check that the number of self-consistent
   iterations is unchanged and that the total energy still matches the baseline
   to the accuracy you care about. A setting that changes the number of
   iterations is changing the calculation, not just its speed, and its timing
   is not comparable.

#. Re-check the winning combination against the baseline in a final set of
   repeats. Gains from separate knobs are often not additive.

#. If the machine has nodes of differing performance, or a shared network, run
   the comparison within a single allocation rather than across separately
   queued jobs, so that the comparison is not confounded by which nodes you
   were given.

Finally, keep the successful settings with the input file. Environment
variables that live only in a submission script are easy to lose, and a
calculation that is 30% slower for no visible reason is hard to diagnose
months later.
