cp ../perfect_density/charges.bin .
dftb+ | tee output
dp_bands band.out band

