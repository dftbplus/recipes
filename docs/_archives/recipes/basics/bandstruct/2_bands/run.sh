#!/bin/bash
set -e
cp ../1_density/charges.bin .
dftb+ > output
dp_bands band.out band >> output

