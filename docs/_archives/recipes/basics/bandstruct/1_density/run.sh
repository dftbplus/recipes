#!/bin/bash
set -e
dftb+ > output
dp_dos band.out dos_total.dat >> output
dp_dos -w dos_ti.1.out dos_ti.s.dat >> output
dp_dos -w dos_ti.2.out dos_ti.p.dat >> output
dp_dos -w dos_ti.3.out dos_ti.d.dat >> output
dp_dos -w dos_o.1.out dos_o.s.dat >> output
dp_dos -w dos_o.2.out dos_o.p.dat >> output

