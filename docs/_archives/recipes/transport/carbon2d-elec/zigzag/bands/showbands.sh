grep "Fermi energy" ../density/detailed.out
plotxy -L --xlabel "K points" --ylabel "Energy [eV]" band_tot.dat & 
