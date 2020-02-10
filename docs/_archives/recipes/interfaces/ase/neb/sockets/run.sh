#!/bin/bash
set -e
python3 ./neb_sockets.py > output
ase gui i2f.traj@-15:
