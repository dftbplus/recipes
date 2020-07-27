grep "Fermi level" detailed.out
plotxy --xlabel "Energy [eV]" --ylabel "DOS" dos.dat pdos.C.1.dat pdos.C.2.dat &
