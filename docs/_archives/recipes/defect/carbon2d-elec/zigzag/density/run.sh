dftb+ | tee output
dp_dos band.out dos.dat
dp_dos -w pdos.C1.out pdos.C1.dat
dp_dos -w pdos.C2.out pdos.C2.dat
dp_dos -w pdos.C3.out pdos.C3.dat
dp_dos -w pdos.H.out pdos.H.dat

