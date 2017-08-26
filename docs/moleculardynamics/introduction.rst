.. highlight:: none

******************
Molecular dynamics
******************

Classical molecular dynamics treats the position and motion of the atoms (the
ion cores) as with Newtonian dynamics. Often this is with modifications to
include the effect of contact with a thermal reservoir (NVT) or additionally
treats the system as though it is in a piston that controls the pressure (NPT).

DFTB+ has several options for basic molecular dynamics, but also can use an
accelerated form of SCC-DFTB which has similar performance properties to non-SCC
for suitable systems.

The basic method adopted for molecular dynamics in DFTB+ is the Velocity Verlet
integrator (Probably more correctly called the Newton-Stormer-Verlet method, as
its been independently discovered several times). This is an implicit symplectic
integrator, which is good at preserving the total energies in the NVE ensemble
for separable Hamilton systems with energy of the form :math:`H(p, q) = T(p) +
U(q)` where :math:`p` and :math:`q` are conjugate momenta and positions of the
atoms.

The form of the integrator is an update for the positions at step :math:`n+1`:

:math:`x_{n+1} = x_n + \tau v_n + \tau^2 F_n / 2 m`,

followed by evaluation of the forces, :math:`F_{n+1}` at that geometry and then
an update for the atomic velocities:

:math:`v_{n+1} = v_n + \tau \left[ F_{n+1} + F_{n} \right] / 2 m`.
      
The user defined parameter :math:`\tau` controls the accuracy of the
integration, typically this should be chosen to be somewhat shorter than the
fastest vibrational period of the system. For many systems, :math:`\tau \approx
1`fs is reasonably stable, but this should of course be tested for the
particular case you are interested in.
