grep "Fermi level" ../density/detailed.out
plotxy -L --xlabel "K points" --ylabel "Energy [eV]" band_tot.dat & 
