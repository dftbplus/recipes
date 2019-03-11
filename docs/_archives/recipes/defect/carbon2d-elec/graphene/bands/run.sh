cp ../latopt/geo_end.gen .
cp ../latopt/charges.bin .
dftb+ | tee output
dp_bands band.out band

