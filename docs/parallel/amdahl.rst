Amdahl's law
============

One of the major issues in the performance of a parallel program is the possible
speed up from using parallelism. Amdahl's law describes the expected gains for
operations which are `partially` parallelised, i.e., where there are
computational parts that must be performed in serial as well others which
operate in parallel. If the fraction of the operations in parallel are `p`, then
compared to serial operation, the expected speed-up is:

.. math::
   
   S(N) = \frac{ 1 }{ 1 - p + \frac{p}{N} }

It can provide a reasonable indication of the performance of a problem on
increasing numbers of processors, depending on the parallelism of the task:

  .. figure:: ../_figures/parallel/amdahl.png
     :height: 40ex
     :align: center
     :alt: Amdahl's law for speed-up.

As you can see, only substantially parallel tasks give substantial speed up. The
limiting speed up is 1 / (1 - p), i.e., the inverse of the serial fraction of
the code.

Amdahl's does not include effects of latencies on parallel performance, but can
be a good guide to the limits of scalability for `fixed` sizes of problem on
parallel machines.
