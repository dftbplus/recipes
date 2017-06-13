#!/bin/bash

(time dftb+)  >& output
cat md.out | awk '/Total MD Energy/ { print $4 }' > energies.dat
velo_autocorr geo_end.xyz xlbomd 0 20000 10000 0.5e-3
